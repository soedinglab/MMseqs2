/**
 * @file test_snappy_diag.c
 * @brief Snappy rejection diagnostics + first-data-page offset resolution.
 *
 * Two things are covered here, both driven by a field report of a 3.6 GB
 * snappy file that failed mid-way through with nothing but
 * CARQUET_ERROR_INVALID_COMPRESSED_DATA to go on:
 *
 *  1. Every distinct reason the Snappy decoder can reject a block must be
 *     reported distinctly, so one failing run on a huge file is enough to
 *     identify the cause.
 *  2. carquet_resolve_data_start_offset() must honour a conforming writer's
 *     data_page_offset and only override it when it is provably wrong.
 */

#include "test_helpers.h"
#include "compression/snappy.h"
#include "reader/reader_internal.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/* ------------------------------------------------------------------ */
/* Helpers                                                             */
/* ------------------------------------------------------------------ */

/* Snappy tag encodings used to hand-build blocks below. */
#define SNAP_TAG_LITERAL(len_minus_1) ((uint8_t)(((len_minus_1) << 2) | 0x00))
#define SNAP_TAG_COPY1(len, off_hi)   ((uint8_t)(0x01 | (((len) - 4) << 2) | ((off_hi) << 5)))
#define SNAP_TAG_COPY2(len)           ((uint8_t)(0x02 | (((len) - 1) << 2)))
#define SNAP_TAG_COPY4(len)           ((uint8_t)(0x03 | (((len) - 1) << 2)))

static size_t put_varint(uint8_t* p, uint32_t v) {
    size_t n = 0;
    while (v >= 0x80) { p[n++] = (uint8_t)(v | 0x80); v >>= 7; }
    p[n++] = (uint8_t)v;
    return n;
}

/* Decompress `src` and assert it is rejected with exactly `want`. */
static void expect_reason(const char* name,
                          const uint8_t* src, size_t src_size,
                          size_t dst_capacity,
                          carquet_snappy_reason_t want) {
    uint8_t* dst = (uint8_t*)malloc(dst_capacity ? dst_capacity : 1);
    assert(dst);
    size_t out = 0;
    carquet_snappy_diag_t diag;
    carquet_status_t st = carquet_snappy_decompress_diag(
        src, src_size, dst, dst_capacity, &out, &diag);

    if (st != CARQUET_ERROR_INVALID_COMPRESSED_DATA) {
        printf("[FAIL] %s: expected rejection, got status %d\n", name, (int)st);
        free(dst);
        exit(1);
    }
    if (diag.reason != want) {
        char buf[192];
        printf("[FAIL] %s: expected reason %d (%s), got %d (%s) -> %s\n",
               name, (int)want, carquet_snappy_reason_string(want),
               (int)diag.reason, carquet_snappy_reason_string(diag.reason),
               carquet_snappy_diag_format(&diag, buf, sizeof(buf)));
        free(dst);
        exit(1);
    }
    /* The formatted line must name the reason, be NUL-terminated, and stay
     * short enough that the reader's own prefix plus this detail still fit in
     * a carquet_error_t message (CARQUET_ERROR_MESSAGE_MAX). */
    char buf[192];
    carquet_snappy_diag_format(&diag, buf, sizeof(buf));
    assert(strstr(buf, "snappy[") != NULL);
    assert(strstr(buf, carquet_snappy_reason_string(want)) != NULL);
    assert(strlen(buf) < sizeof(buf));
    assert(strlen(buf) <= CARQUET_ERROR_MESSAGE_MAX - 100);

    free(dst);
    TEST_PASS(name);
}

/* ------------------------------------------------------------------ */
/* 1. Distinct rejection reasons                                       */
/* ------------------------------------------------------------------ */

static int test_reason_length_varint(void) {
    /* Five continuation bytes: the length varint cannot terminate. */
    const uint8_t src[] = {0x80, 0x80, 0x80, 0x80, 0x80};
    expect_reason("snappy_diag/length_varint", src, sizeof(src), 64,
                  CARQUET_SNAPPY_ERR_LENGTH_VARINT);
    return 0;
}

static int test_reason_output_too_small(void) {
    uint8_t src[16];
    size_t n = put_varint(src, 1000);
    src[n++] = SNAP_TAG_LITERAL(0);
    src[n++] = 'x';
    expect_reason("snappy_diag/output_too_small", src, n, 8,
                  CARQUET_SNAPPY_ERR_OUTPUT_TOO_SMALL);
    return 0;
}

static int test_reason_literal_len_trunc(void) {
    uint8_t src[8];
    size_t n = put_varint(src, 100);
    src[n++] = SNAP_TAG_LITERAL(62);  /* 62 -> 3 extra length bytes follow */
    src[n++] = 0x01;             /* only one of them present */
    expect_reason("snappy_diag/literal_len_trunc", src, n, 128,
                  CARQUET_SNAPPY_ERR_LITERAL_LEN_TRUNC);
    return 0;
}

static int test_reason_literal_input(void) {
    /* Declares a 40-byte literal but supplies 4 bytes. Kept above the
     * 16-byte fast path so the bounds check is the one that fires. */
    uint8_t src[16];
    size_t n = put_varint(src, 40);
    src[n++] = SNAP_TAG_LITERAL(39);
    src[n++] = 'a'; src[n++] = 'b'; src[n++] = 'c'; src[n++] = 'd';
    expect_reason("snappy_diag/literal_input", src, n, 128,
                  CARQUET_SNAPPY_ERR_LITERAL_INPUT);
    return 0;
}

static int test_reason_literal_output(void) {
    /* Header says 4 bytes; the literal claims 40 and the input has them. */
    uint8_t src[64];
    size_t n = put_varint(src, 4);
    src[n++] = SNAP_TAG_LITERAL(39);
    for (int i = 0; i < 40; i++) src[n++] = (uint8_t)('a' + (i % 26));
    expect_reason("snappy_diag/literal_output", src, n, 128,
                  CARQUET_SNAPPY_ERR_LITERAL_OUTPUT);
    return 0;
}

static int test_reason_copy1_trunc(void) {
    uint8_t src[64];
    size_t n = put_varint(src, 100);
    src[n++] = SNAP_TAG_LITERAL(19);
    for (int i = 0; i < 20; i++) src[n++] = 'z';
    src[n++] = SNAP_TAG_COPY1(4, 0);  /* trailer byte missing */
    expect_reason("snappy_diag/copy1_trunc", src, n, 256,
                  CARQUET_SNAPPY_ERR_COPY1_TRUNC);
    return 0;
}

static int test_reason_copy2_trunc(void) {
    uint8_t src[64];
    size_t n = put_varint(src, 100);
    src[n++] = SNAP_TAG_LITERAL(19);
    for (int i = 0; i < 20; i++) src[n++] = 'z';
    src[n++] = SNAP_TAG_COPY2(8);
    src[n++] = 0x01;             /* only one of the two trailer bytes */
    expect_reason("snappy_diag/copy2_trunc", src, n, 256,
                  CARQUET_SNAPPY_ERR_COPY2_TRUNC);
    return 0;
}

static int test_reason_copy4_trunc(void) {
    uint8_t src[64];
    size_t n = put_varint(src, 100);
    src[n++] = SNAP_TAG_LITERAL(19);
    for (int i = 0; i < 20; i++) src[n++] = 'z';
    src[n++] = SNAP_TAG_COPY4(8);
    src[n++] = 0x01; src[n++] = 0x00;  /* two of the four trailer bytes */
    expect_reason("snappy_diag/copy4_trunc", src, n, 256,
                  CARQUET_SNAPPY_ERR_COPY4_TRUNC);
    return 0;
}

static int test_reason_copy_offset_zero(void) {
    uint8_t src[64];
    size_t n = put_varint(src, 100);
    src[n++] = SNAP_TAG_LITERAL(19);
    for (int i = 0; i < 20; i++) src[n++] = 'z';
    src[n++] = SNAP_TAG_COPY2(8);
    src[n++] = 0x00; src[n++] = 0x00;  /* offset 0 */
    expect_reason("snappy_diag/copy_offset_zero", src, n, 256,
                  CARQUET_SNAPPY_ERR_COPY_OFFSET_ZERO);
    return 0;
}

static int test_reason_copy_offset_range(void) {
    uint8_t src[64];
    size_t n = put_varint(src, 100);
    src[n++] = SNAP_TAG_LITERAL(19);
    for (int i = 0; i < 20; i++) src[n++] = 'z';
    src[n++] = SNAP_TAG_COPY2(8);
    src[n++] = 0x64; src[n++] = 0x00;  /* offset 100 > 20 bytes produced */
    expect_reason("snappy_diag/copy_offset_range", src, n, 256,
                  CARQUET_SNAPPY_ERR_COPY_OFFSET_RANGE);
    return 0;
}

static int test_reason_copy_output(void) {
    /* 20 literal bytes then a 64-byte copy, but only 24 bytes declared. */
    uint8_t src[64];
    size_t n = put_varint(src, 24);
    src[n++] = SNAP_TAG_LITERAL(19);
    for (int i = 0; i < 20; i++) src[n++] = 'z';
    src[n++] = SNAP_TAG_COPY2(64);
    src[n++] = 0x04; src[n++] = 0x00;
    expect_reason("snappy_diag/copy_output", src, n, 256,
                  CARQUET_SNAPPY_ERR_COPY_OUTPUT);
    return 0;
}

static int test_reason_short_output(void) {
    /* Declares 100 bytes, delivers 20 and then runs out of input. */
    uint8_t src[64];
    size_t n = put_varint(src, 100);
    src[n++] = SNAP_TAG_LITERAL(19);
    for (int i = 0; i < 20; i++) src[n++] = 'z';
    expect_reason("snappy_diag/short_output", src, n, 256,
                  CARQUET_SNAPPY_ERR_SHORT_OUTPUT);
    return 0;
}

static int test_reason_trailing_input(void) {
    /* Output is complete after 20 bytes but three stray bytes remain. */
    uint8_t src[64];
    size_t n = put_varint(src, 20);
    src[n++] = SNAP_TAG_LITERAL(19);
    for (int i = 0; i < 20; i++) src[n++] = 'z';
    src[n++] = 0x00; src[n++] = 0x00; src[n++] = 0x00;
    expect_reason("snappy_diag/trailing_input", src, n, 256,
                  CARQUET_SNAPPY_ERR_TRAILING_INPUT);
    return 0;
}

/* Every reason must map to its own string, and no two may collide. */
static int test_reason_strings_unique(void) {
    const carquet_snappy_reason_t all[] = {
        CARQUET_SNAPPY_ERR_LENGTH_VARINT, CARQUET_SNAPPY_ERR_OUTPUT_TOO_SMALL,
        CARQUET_SNAPPY_ERR_LITERAL_LEN_TRUNC, CARQUET_SNAPPY_ERR_LITERAL_INPUT,
        CARQUET_SNAPPY_ERR_LITERAL_OUTPUT, CARQUET_SNAPPY_ERR_COPY1_TRUNC,
        CARQUET_SNAPPY_ERR_COPY2_TRUNC, CARQUET_SNAPPY_ERR_COPY4_TRUNC,
        CARQUET_SNAPPY_ERR_COPY_OFFSET_ZERO, CARQUET_SNAPPY_ERR_COPY_OFFSET_RANGE,
        CARQUET_SNAPPY_ERR_COPY_OUTPUT, CARQUET_SNAPPY_ERR_SHORT_OUTPUT,
        CARQUET_SNAPPY_ERR_TRAILING_INPUT,
    };
    const size_t n = sizeof(all) / sizeof(all[0]);
    for (size_t i = 0; i < n; i++) {
        const char* a = carquet_snappy_reason_string(all[i]);
        assert(a && *a);
        assert(strcmp(a, "unknown") != 0);
        for (size_t j = i + 1; j < n; j++) {
            if (strcmp(a, carquet_snappy_reason_string(all[j])) == 0) {
                TEST_FAIL("snappy_diag/reason_strings_unique", "duplicate text");
            }
        }
    }
    TEST_PASS("snappy_diag/reason_strings_unique");
    return 0;
}

/* A valid block must decode cleanly and leave the diagnostic untouched. */
static int test_valid_block_reports_ok(void) {
    uint8_t plain[4096];
    for (size_t i = 0; i < sizeof(plain); i++) {
        plain[i] = (uint8_t)((i / 7) ^ (i % 13));
    }
    size_t bound = carquet_snappy_compress_bound(sizeof(plain));
    uint8_t* comp = (uint8_t*)malloc(bound);
    assert(comp);
    size_t comp_size = 0;
    assert(carquet_snappy_compress(plain, sizeof(plain), comp, bound, &comp_size)
           == CARQUET_OK);

    uint8_t out[sizeof(plain)];
    size_t out_size = 0;
    carquet_snappy_diag_t diag;
    assert(carquet_snappy_decompress_diag(comp, comp_size, out, sizeof(out),
                                          &out_size, &diag) == CARQUET_OK);
    assert(out_size == sizeof(plain));
    assert(memcmp(out, plain, sizeof(plain)) == 0);
    assert(diag.reason == CARQUET_SNAPPY_OK);

    /* The plain entry point must stay behaviourally identical. */
    size_t out_size2 = 0;
    assert(carquet_snappy_decompress(comp, comp_size, out, sizeof(out), &out_size2)
           == CARQUET_OK);
    assert(out_size2 == sizeof(plain));

    free(comp);
    TEST_PASS("snappy_diag/valid_block");
    return 0;
}

/* ------------------------------------------------------------------ */
/* 2. First-data-page offset resolution                                */
/* ------------------------------------------------------------------ */

static int test_data_start_offset_resolution(void) {
    /* Conforming writer: the dictionary page is immediately followed by the
     * first data page, and both agree. */
    assert(carquet_resolve_data_start_offset(1000, 1000) == 1000);

    /* Conforming writer with something in between (index page, alignment
     * padding): the writer's offset must win. Recomputing here is what made
     * carquet read a data page out of the middle of another page. */
    assert(carquet_resolve_data_start_offset(1064, 1000) == 1064);
    assert(carquet_resolve_data_start_offset(9000, 1000) == 9000);

    /* Broken writer (DuckDB): data_page_offset points at the dictionary page,
     * i.e. before its end. Fall back to the computed end of the dictionary. */
    assert(carquet_resolve_data_start_offset(900, 1000) == 1000);
    assert(carquet_resolve_data_start_offset(0, 1000) == 1000);
    assert(carquet_resolve_data_start_offset(999, 1000) == 1000);

    /* Offsets past 2 GiB must be handled as 64-bit throughout. */
    const int64_t big = (int64_t)3 * 1024 * 1024 * 1024;
    assert(carquet_resolve_data_start_offset(big + 4096, big) == big + 4096);
    assert(carquet_resolve_data_start_offset(big - 10, big) == big);

    TEST_PASS("snappy_diag/data_start_offset_resolution");
    return 0;
}

/* ------------------------------------------------------------------ */
/* 3. The diagnostic must survive out to the batch-reader API           */
/* ------------------------------------------------------------------ */

#define BR_ROWS 20000
#define BR_COLS 3

/* Write a 3-column snappy file; return the first data page offset of column 0
 * through @p data_page_offset. */
static int write_snappy_file(const char* path, int64_t* data_page_offset) {
    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema) return -1;
    for (int c = 0; c < BR_COLS; c++) {
        char name[16];
        snprintf(name, sizeof(name), "s%d", c);
        carquet_schema_add_column(schema, name, CARQUET_PHYSICAL_BYTE_ARRAY, NULL,
                                  CARQUET_REPETITION_REQUIRED, 0, 0);
    }

    carquet_writer_options_t wo;
    carquet_writer_options_init(&wo);
    wo.compression = CARQUET_COMPRESSION_SNAPPY;
    /* Several row groups so the multi-row-group pipeline path is exercised. */
    wo.row_group_size = BR_ROWS / 4;
    carquet_writer_t* w = carquet_writer_create(path, schema, &wo, &err);
    if (!w) { carquet_schema_free(schema); return -1; }

    carquet_byte_array_t* vals = (carquet_byte_array_t*)malloc(sizeof(*vals) * BR_ROWS);
    char* store = (char*)malloc((size_t)BR_ROWS * 64);
    assert(vals && store);
    for (int i = 0; i < BR_ROWS; i++) {
        char* p = store + (size_t)i * 64;
        int n = snprintf(p, 64, "row-%08d-payload-abcdefghijklmnopqrstuvwxyz", i);
        vals[i].data = (uint8_t*)p;
        vals[i].length = (uint32_t)n;
    }
    for (int c = 0; c < BR_COLS; c++) {
        carquet_writer_write_batch(w, c, vals, BR_ROWS, NULL, NULL);
    }
    carquet_writer_close(w);
    free(vals);
    free(store);
    carquet_schema_free(schema);

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    if (!r) return -1;
    carquet_column_chunk_metadata_t m;
    if (carquet_reader_column_chunk_metadata(r, 0, 0, &m) != CARQUET_OK) {
        carquet_reader_close(r);
        return -1;
    }
    *data_page_offset = m.has_dictionary_page && m.dictionary_page_offset < m.data_page_offset
                      ? m.dictionary_page_offset : m.data_page_offset;
    carquet_reader_close(r);
    return 0;
}

/* Overwrite @p len bytes at @p offset so the page payload stops being a valid
 * Snappy block without disturbing the page header. */
static void corrupt_at(const char* path, int64_t offset, int len) {
    FILE* f = fopen(path, "r+b");
    assert(f);
    assert(fseek(f, (long)offset, SEEK_SET) == 0);
    for (int i = 0; i < len; i++) fputc(0xA5, f);
    fclose(f);
}

static carquet_batch_reader_t* open_batch_reader(const char* path,
                                                 bool use_mmap,
                                                 carquet_reader_t** out_reader) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_options_t ro;
    carquet_reader_options_init(&ro);
    /* Reach the codec: with CRC verification on, the page checksum rejects the
     * corruption first and we never exercise the decompression path. */
    ro.verify_checksums = false;
    ro.use_mmap = use_mmap;
    carquet_reader_t* r = carquet_reader_open(path, &ro, &err);
    if (!r) return NULL;
    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    cfg.batch_size = 4096;
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    *out_reader = r;
    return br;
}

static int test_batch_reader_last_error(void) {
    /* NULL reader is accepted and reports nothing. */
    assert(carquet_batch_reader_last_error(NULL) == NULL);

    const char* path = "test_snappy_diag_br.parquet";
    int64_t page_offset = 0;
    if (write_snappy_file(path, &page_offset) != 0) {
        TEST_FAIL("snappy_diag/batch_reader_last_error", "failed to write fixture");
    }

    /* --- Clean file: every successful next() leaves the error unset, and so
     *     does the END_OF_DATA terminator. --- */
    carquet_reader_t* r = NULL;
    carquet_batch_reader_t* br = open_batch_reader(path, false, &r);
    assert(br);
    carquet_row_batch_t* b = NULL;
    carquet_status_t st;
    int64_t total = 0;
    while ((st = carquet_batch_reader_next(br, &b)) == CARQUET_OK && b) {
        assert(carquet_batch_reader_last_error(br) == NULL);
        total += carquet_row_batch_num_rows(b);
        carquet_row_batch_free(b);
        b = NULL;
    }
    assert(total == BR_ROWS);
    assert(st == CARQUET_ERROR_END_OF_DATA || st == CARQUET_OK);
    assert(carquet_batch_reader_last_error(br) == NULL);
    carquet_batch_reader_free(br);
    carquet_reader_close(r);

    /* --- Corrupted payload: next() returns a status, and last_error carries
     *     the page offset plus the exact Snappy check that fired. --- */
    corrupt_at(path, page_offset + 200, 64);

    /* Both I/O modes: fread and mmap (which also drives the parallel
     * multi-column prefetch, so the per-column error slots are exercised). */
    for (int mmap_mode = 0; mmap_mode < 2; mmap_mode++) {
        r = NULL;
        br = open_batch_reader(path, mmap_mode != 0, &r);
        assert(br);
        b = NULL;
        /* The corrupted page is not necessarily in the first batch. */
        int failed = 0;
        while ((st = carquet_batch_reader_next(br, &b)) == CARQUET_OK && b) {
            carquet_row_batch_free(b);
            b = NULL;
        }
        failed = (st != CARQUET_OK && st != CARQUET_ERROR_END_OF_DATA);
        if (!failed) {
            carquet_batch_reader_free(br);
            carquet_reader_close(r);
            remove(path);
            TEST_FAIL("snappy_diag/batch_reader_last_error",
                      "corrupted page still read cleanly");
        }

        const carquet_error_t* e = carquet_batch_reader_last_error(br);
        if (!e) {
            carquet_batch_reader_free(br);
            carquet_reader_close(r);
            remove(path);
            TEST_FAIL("snappy_diag/batch_reader_last_error",
                      "failure reported no detail");
        }
        /* The underlying codec status, not the coarser status next() returns. */
        assert(e->code == CARQUET_ERROR_INVALID_COMPRESSED_DATA);
        assert(strstr(e->message, "snappy[") != NULL);
        assert(strstr(e->message, "file offset") != NULL);
        /* Nothing was truncated away by CARQUET_ERROR_MESSAGE_MAX. */
        assert(e->message[0] != '\0');
        assert(strlen(e->message) < CARQUET_ERROR_MESSAGE_MAX - 1);

        carquet_batch_reader_free(br);
        carquet_reader_close(r);
    }
    remove(path);

    TEST_PASS("snappy_diag/batch_reader_last_error");
    return 0;
}

int main(void) {
    printf("=== Snappy diagnostics and page offset resolution ===\n");

    if (test_reason_length_varint()) return 1;
    if (test_reason_output_too_small()) return 1;
    if (test_reason_literal_len_trunc()) return 1;
    if (test_reason_literal_input()) return 1;
    if (test_reason_literal_output()) return 1;
    if (test_reason_copy1_trunc()) return 1;
    if (test_reason_copy2_trunc()) return 1;
    if (test_reason_copy4_trunc()) return 1;
    if (test_reason_copy_offset_zero()) return 1;
    if (test_reason_copy_offset_range()) return 1;
    if (test_reason_copy_output()) return 1;
    if (test_reason_short_output()) return 1;
    if (test_reason_trailing_input()) return 1;
    if (test_reason_strings_unique()) return 1;
    if (test_valid_block_reports_ok()) return 1;
    if (test_data_start_offset_resolution()) return 1;
    if (test_batch_reader_last_error()) return 1;

    printf("=== All snappy diagnostic tests passed ===\n");
    return 0;
}
