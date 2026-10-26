/**
 * @file test_parallel_write.c
 * @brief Page-granular parallel encode + background I/O in the writer.
 *
 * The row-group finalize schedules one task per page for flat fixed-stride
 * columns (src/writer/page_tasks.h) and hands finished row groups to an I/O
 * thread (carquet_writer_options_t.async_io). Neither may change a single
 * output byte: this test writes the same columns through the serial path
 * (a one-column file never finalizes in parallel) and through the parallel
 * path, then compares the column chunk bytes, and writes the same file with
 * async_io on and off and compares the whole files. It also drives several
 * row groups through the buffer recycling of the I/O thread and reads every
 * value back.
 */

#include <carquet/carquet.h>
#include "test_helpers.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define ROWS 300000   /* > 2 pages of 8-byte values at the 1 MB default */

static uint32_t lcg_state;
static uint32_t lcg(void) {
    lcg_state = lcg_state * 1103515245u + 12345u;
    return (lcg_state >> 8) & 0xFFFFFFu;
}

typedef struct {
    int64_t* ids;
    double* values;
    int32_t* cats;
} data_t;

static void data_init(data_t* d) {
    d->ids = malloc(ROWS * sizeof(int64_t));
    d->values = malloc(ROWS * sizeof(double));
    d->cats = malloc(ROWS * sizeof(int32_t));
    lcg_state = 7;
    for (int i = 0; i < ROWS; i++) {
        d->ids[i] = 1000000 + (int64_t)lcg();
        d->values[i] = (double)lcg() / 1024.0 - 3000.0;
        d->cats[i] = (int32_t)(lcg() % 97);
    }
}

static void data_free(data_t* d) {
    free(d->ids); free(d->values); free(d->cats);
}

typedef struct {
    carquet_compression_t codec;
    int64_t write_batch_size;
    int64_t max_rows_per_page;
    int data_page_version;
    bool async_io;
    int64_t rg_rows;   /* 0 = one row group */
} cfg_t;

static void apply_cfg(carquet_writer_options_t* o, const cfg_t* c) {
    carquet_writer_options_init(o);
    o->compression = c->codec;
    o->compression_level = c->codec == CARQUET_COMPRESSION_ZSTD ? 1 : 0;
    o->write_batch_size = c->write_batch_size;
    o->max_rows_per_page = c->max_rows_per_page;
    o->data_page_version = c->data_page_version;
    o->async_io = c->async_io;
}

/* Write `ncols` of the three columns (in order id, value, cat). */
static int write_file(const char* path, const data_t* d, int ncols, const cfg_t* c) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* s = carquet_schema_create(&err);
    if (!s) return 0;
    static const char* names[3] = {"id", "value", "cat"};
    static const carquet_physical_type_t types[3] = {
        CARQUET_PHYSICAL_INT64, CARQUET_PHYSICAL_DOUBLE, CARQUET_PHYSICAL_INT32};
    for (int i = 0; i < ncols; i++) {
        if (carquet_schema_add_column(s, names[i], types[i], NULL,
                CARQUET_REPETITION_REQUIRED, 0, 0) != CARQUET_OK) {
            carquet_schema_free(s); return 0;
        }
    }
    carquet_writer_options_t o; apply_cfg(&o, c);
    carquet_writer_t* w = carquet_writer_create(path, s, &o, &err);
    if (!w) { carquet_schema_free(s); return 0; }
    int ok = 1;
    int64_t step = c->rg_rows > 0 ? c->rg_rows : ROWS;
    for (int64_t off = 0; off < ROWS && ok; off += step) {
        int64_t n = ROWS - off < step ? ROWS - off : step;
        const void* cols[3] = {d->ids + off, d->values + off, d->cats + off};
        for (int i = 0; i < ncols && ok; i++) {
            if (carquet_writer_write_batch(w, i, cols[i], n, NULL, NULL) != CARQUET_OK) ok = 0;
        }
        if (ok && off + n < ROWS && carquet_writer_new_row_group(w) != CARQUET_OK) ok = 0;
    }
    if (carquet_writer_close(w) != CARQUET_OK) ok = 0;
    carquet_schema_free(s);
    return ok;
}

static uint8_t* read_range(const char* path, int64_t off, int64_t len) {
    FILE* f = fopen(path, "rb");
    if (!f) return NULL;
    uint8_t* buf = malloc((size_t)len);
    if (!buf || fseek(f, (long)off, SEEK_SET) != 0 ||
        fread(buf, 1, (size_t)len, f) != (size_t)len) {
        free(buf); buf = NULL;
    }
    fclose(f);
    return buf;
}

static long file_size(const char* path) {
    FILE* f = fopen(path, "rb");
    if (!f) return -1;
    fseek(f, 0, SEEK_END);
    long n = ftell(f);
    fclose(f);
    return n;
}

/* Compare column `col` of the 3-column file against column 0 of a file that
 * holds only that column: same chunk size, same bytes. */
static int chunk_matches(const char* multi, const char* single, int col, const char* tag) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* rm = carquet_reader_open(multi, NULL, &err);
    carquet_reader_t* rs = carquet_reader_open(single, NULL, &err);
    if (!rm || !rs) {
        if (rm) carquet_reader_close(rm);
        if (rs) carquet_reader_close(rs);
        printf("  [%s] open failed\n", tag); return 0;
    }
    int ok = 0;
    carquet_column_chunk_metadata_t mm, ms;
    if (carquet_reader_column_chunk_metadata(rm, 0, col, &mm) == CARQUET_OK &&
        carquet_reader_column_chunk_metadata(rs, 0, 0, &ms) == CARQUET_OK) {
        if (mm.total_compressed_size != ms.total_compressed_size ||
            mm.num_values != ms.num_values || mm.num_values != ROWS) {
            printf("  [%s] col %d size/count differ: %lld/%lld vs %lld/%lld\n", tag, col,
                   (long long)mm.total_compressed_size, (long long)mm.num_values,
                   (long long)ms.total_compressed_size, (long long)ms.num_values);
        } else {
            uint8_t* a = read_range(multi, mm.data_page_offset, mm.total_compressed_size);
            uint8_t* b = read_range(single, ms.data_page_offset, ms.total_compressed_size);
            ok = a && b && memcmp(a, b, (size_t)mm.total_compressed_size) == 0;
            if (!ok) printf("  [%s] col %d chunk bytes differ\n", tag, col);
            free(a); free(b);
        }
    }
    carquet_reader_close(rm);
    carquet_reader_close(rs);
    return ok;
}

static int values_match(const char* path, const data_t* d, int ncols, const char* tag) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_options_t ro; carquet_reader_options_init(&ro);
    ro.use_mmap = true;
    carquet_reader_t* r = carquet_reader_open(path, &ro, &err);
    if (!r) { printf("  [%s] open failed\n", tag); return 0; }
    carquet_batch_reader_config_t c; carquet_batch_reader_config_init(&c);
    c.batch_size = 65536;
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &c, &err);
    int ok = br != NULL;
    int64_t row = 0;
    carquet_row_batch_t* b = NULL;
    while (ok && carquet_batch_reader_next(br, &b) == CARQUET_OK && b) {
        int64_t n = carquet_row_batch_num_rows(b);
        const void* p; const uint8_t* nulls; int64_t cnt;
        if (ncols > 0 && carquet_row_batch_column(b, 0, &p, &nulls, &cnt) == CARQUET_OK) {
            if (cnt != n || memcmp(p, d->ids + row, (size_t)n * 8) != 0) ok = 0;
        }
        if (ncols > 1 && carquet_row_batch_column(b, 1, &p, &nulls, &cnt) == CARQUET_OK) {
            if (cnt != n || memcmp(p, d->values + row, (size_t)n * 8) != 0) ok = 0;
        }
        if (ncols > 2 && carquet_row_batch_column(b, 2, &p, &nulls, &cnt) == CARQUET_OK) {
            if (cnt != n || memcmp(p, d->cats + row, (size_t)n * 4) != 0) ok = 0;
        }
        row += n;
        carquet_row_batch_free(b); b = NULL;
    }
    if (row != ROWS) ok = 0;
    if (!ok) printf("  [%s] value mismatch (rows=%lld)\n", tag, (long long)row);
    if (br) carquet_batch_reader_free(br);
    carquet_reader_close(r);
    return ok;
}

static int test_page_parallel_identical(const cfg_t* c, const char* tag) {
    data_t d; data_init(&d);
    char multi[512], single[3][512];
    carquet_test_temp_path(multi, sizeof(multi), "ppw_multi");
    carquet_test_temp_path(single[0], sizeof(single[0]), "ppw_s0");
    carquet_test_temp_path(single[1], sizeof(single[1]), "ppw_s1");
    carquet_test_temp_path(single[2], sizeof(single[2]), "ppw_s2");

    int ok = write_file(multi, &d, 3, c);
    /* Column i alone: the one-column row group finalizes serially. */
    for (int i = 0; ok && i < 3; i++) {
        /* Build the one-column file by hand to keep the physical type. */
        carquet_error_t err = CARQUET_ERROR_INIT;
        carquet_schema_t* s = carquet_schema_create(&err);
        static const carquet_physical_type_t types[3] = {
            CARQUET_PHYSICAL_INT64, CARQUET_PHYSICAL_DOUBLE, CARQUET_PHYSICAL_INT32};
        ok = s && carquet_schema_add_column(s, "c", types[i], NULL,
                CARQUET_REPETITION_REQUIRED, 0, 0) == CARQUET_OK;
        carquet_writer_options_t o; apply_cfg(&o, c);
        carquet_writer_t* w = ok ? carquet_writer_create(single[i], s, &o, &err) : NULL;
        ok = w != NULL;
        const void* src = i == 0 ? (const void*)d.ids : i == 1 ? (const void*)d.values : (const void*)d.cats;
        if (ok && carquet_writer_write_batch(w, 0, src, ROWS, NULL, NULL) != CARQUET_OK) ok = 0;
        if (w && carquet_writer_close(w) != CARQUET_OK) ok = 0;
        if (s) carquet_schema_free(s);
    }
    for (int i = 0; ok && i < 3; i++) ok = chunk_matches(multi, single[i], i, tag);
    if (ok) ok = values_match(multi, &d, 3, tag);

    carquet_test_cleanup(multi);
    for (int i = 0; i < 3; i++) carquet_test_cleanup(single[i]);
    data_free(&d);
    if (!ok) TEST_FAIL(tag, "parallel page output differs from serial output");
    TEST_PASS(tag);
    return 0;
}

static int test_async_io_identical(carquet_compression_t codec, const char* tag) {
    data_t d; data_init(&d);
    char a[512], b[512];
    carquet_test_temp_path(a, sizeof(a), "aio_on");
    carquet_test_temp_path(b, sizeof(b), "aio_off");
    cfg_t on = {codec, 0, 0, 1, true, 70000};
    cfg_t off = on; off.async_io = false;
    int ok = write_file(a, &d, 3, &on) && write_file(b, &d, 3, &off);
    if (ok) {
        long la = file_size(a), lb = file_size(b);
        ok = la > 0 && la == lb;
        if (ok) {
            uint8_t* pa = read_range(a, 0, la);
            uint8_t* pb = read_range(b, 0, lb);
            ok = pa && pb && memcmp(pa, pb, (size_t)la) == 0;
            free(pa); free(pb);
        }
        if (!ok) printf("  [%s] files differ (%ld vs %ld bytes)\n", tag, la, lb);
    }
    if (ok) ok = values_match(a, &d, 3, tag);
    carquet_test_cleanup(a); carquet_test_cleanup(b);
    data_free(&d);
    if (!ok) TEST_FAIL(tag, "async_io changed the output");
    TEST_PASS(tag);
    return 0;
}

/* Several row groups through the buffer writer: the I/O thread writes into
 * the temp stream while later groups are encoded, and get_buffer must see
 * the complete file. */
static int test_async_buffer_writer(void) {
    data_t d; data_init(&d);
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* s = carquet_schema_create(&err);
    int ok = s &&
        carquet_schema_add_column(s, "id", CARQUET_PHYSICAL_INT64, NULL,
            CARQUET_REPETITION_REQUIRED, 0, 0) == CARQUET_OK &&
        carquet_schema_add_column(s, "value", CARQUET_PHYSICAL_DOUBLE, NULL,
            CARQUET_REPETITION_REQUIRED, 0, 0) == CARQUET_OK;
    cfg_t c = {CARQUET_COMPRESSION_SNAPPY, 0, 0, 1, true, 40000};
    carquet_writer_options_t o; apply_cfg(&o, &c);
    carquet_writer_t* w = ok ? carquet_writer_create_buffer(s, &o, &err) : NULL;
    ok = w != NULL;
    for (int64_t off = 0; ok && off < ROWS; off += c.rg_rows) {
        int64_t n = ROWS - off < c.rg_rows ? ROWS - off : c.rg_rows;
        ok = carquet_writer_write_batch(w, 0, d.ids + off, n, NULL, NULL) == CARQUET_OK &&
             carquet_writer_write_batch(w, 1, d.values + off, n, NULL, NULL) == CARQUET_OK;
        if (ok && off + n < ROWS) ok = carquet_writer_new_row_group(w) == CARQUET_OK;
    }
    if (ok) ok = carquet_writer_close(w) == CARQUET_OK;
    void* buf = NULL; size_t size = 0;
    if (ok) ok = carquet_writer_get_buffer(w, &buf, &size) == CARQUET_OK && buf && size > 0;
    else if (w) carquet_writer_abort(w);
    if (s) carquet_schema_free(s);

    if (ok) {
        char path[512]; carquet_test_temp_path(path, sizeof(path), "aio_buf");
        FILE* f = fopen(path, "wb");
        ok = f && fwrite(buf, 1, size, f) == size;
        if (f) fclose(f);
        if (ok) ok = values_match(path, &d, 2, "async_buffer_writer");
        carquet_test_cleanup(path);
    }
    free(buf);
    data_free(&d);
    if (!ok) TEST_FAIL("async_buffer_writer", "multi-row-group buffer write failed");
    TEST_PASS("async_buffer_writer");
    return 0;
}

/* A nullable column keeps the serial per-column path while its neighbours
 * take the page-parallel one; both must land correctly in the same file. */
static int test_mixed_nullable(void) {
    data_t d; data_init(&d);
    carquet_error_t err = CARQUET_ERROR_INIT;
    char path[512]; carquet_test_temp_path(path, sizeof(path), "ppw_mixed");
    carquet_schema_t* s = carquet_schema_create(&err);
    int ok = s &&
        carquet_schema_add_column(s, "id", CARQUET_PHYSICAL_INT64, NULL,
            CARQUET_REPETITION_REQUIRED, 0, 0) == CARQUET_OK &&
        carquet_schema_add_column(s, "opt", CARQUET_PHYSICAL_INT32, NULL,
            CARQUET_REPETITION_OPTIONAL, 0, 0) == CARQUET_OK &&
        carquet_schema_add_column(s, "value", CARQUET_PHYSICAL_DOUBLE, NULL,
            CARQUET_REPETITION_REQUIRED, 0, 0) == CARQUET_OK;
    int16_t* defs = malloc(ROWS * sizeof(int16_t));
    int32_t* packed = malloc(ROWS * sizeof(int32_t));
    int64_t np = 0;
    for (int i = 0; i < ROWS; i++) {
        defs[i] = (i % 5 == 0) ? 0 : 1;
        if (defs[i]) packed[np++] = d.cats[i];
    }
    cfg_t c = {CARQUET_COMPRESSION_LZ4_RAW, 0, 0, 1, true, 0};
    carquet_writer_options_t o; apply_cfg(&o, &c);
    carquet_writer_t* w = ok ? carquet_writer_create(path, s, &o, &err) : NULL;
    ok = w != NULL;
    if (ok) ok = carquet_writer_write_batch(w, 0, d.ids, ROWS, NULL, NULL) == CARQUET_OK &&
                 carquet_writer_write_batch(w, 1, packed, ROWS, defs, NULL) == CARQUET_OK &&
                 carquet_writer_write_batch(w, 2, d.values, ROWS, NULL, NULL) == CARQUET_OK;
    if (w && carquet_writer_close(w) != CARQUET_OK) ok = 0;
    if (s) carquet_schema_free(s);

    if (ok) {
        carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
        ok = r != NULL;
        if (ok) {
            int64_t* ids = malloc(ROWS * sizeof(int64_t));
            double* vals = malloc(ROWS * sizeof(double));
            int32_t* got = malloc(ROWS * sizeof(int32_t));
            int16_t* gdef = malloc(ROWS * sizeof(int16_t));
            carquet_column_reader_t* c0 = carquet_reader_get_column(r, 0, 0, &err);
            carquet_column_reader_t* c1 = carquet_reader_get_column(r, 0, 1, &err);
            carquet_column_reader_t* c2 = carquet_reader_get_column(r, 0, 2, &err);
            ok = c0 && c1 && c2 &&
                 carquet_column_read_batch(c0, ids, ROWS, NULL, NULL) == ROWS &&
                 carquet_column_read_batch(c1, got, ROWS, gdef, NULL) == ROWS &&
                 carquet_column_read_batch(c2, vals, ROWS, NULL, NULL) == ROWS &&
                 memcmp(ids, d.ids, ROWS * 8) == 0 &&
                 memcmp(vals, d.values, ROWS * 8) == 0 &&
                 memcmp(gdef, defs, ROWS * sizeof(int16_t)) == 0 &&
                 memcmp(got, packed, (size_t)np * sizeof(int32_t)) == 0;
            if (c0) carquet_column_reader_free(c0);
            if (c1) carquet_column_reader_free(c1);
            if (c2) carquet_column_reader_free(c2);
            free(ids); free(vals); free(got); free(gdef);
            carquet_reader_close(r);
        }
    }
    carquet_test_cleanup(path);
    free(defs); free(packed);
    data_free(&d);
    if (!ok) TEST_FAIL("mixed_nullable", "roundtrip failed");
    TEST_PASS("mixed_nullable");
    return 0;
}

/* A background write failure must surface from close() rather than vanish.
 * Only runnable where a stream that rejects writes exists (/dev/full). */
static int test_async_write_error(void) {
#ifdef _WIN32
    /* Windows has no always-full stream. Do not even try the fopen: under
     * Wine "/dev/full" resolves through the Z: drive to the real Unix device
     * and the UCRT terminates the process inside fopen() with exit status 0,
     * which would silently cut the run short and hide every later test. */
    TEST_PASS("async_write_error (skipped: no /dev/full on Windows)");
    return 0;
#else
    FILE* full = fopen("/dev/full", "wb");
    if (!full) { TEST_PASS("async_write_error (skipped: no /dev/full)"); return 0; }
    data_t d; data_init(&d);
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* s = carquet_schema_create(&err);
    int ok = s &&
        carquet_schema_add_column(s, "id", CARQUET_PHYSICAL_INT64, NULL,
            CARQUET_REPETITION_REQUIRED, 0, 0) == CARQUET_OK &&
        carquet_schema_add_column(s, "value", CARQUET_PHYSICAL_DOUBLE, NULL,
            CARQUET_REPETITION_REQUIRED, 0, 0) == CARQUET_OK;
    carquet_writer_options_t o; carquet_writer_options_init(&o);
    carquet_writer_t* w = ok ? carquet_writer_create_file(full, s, &o, &err) : NULL;
    ok = w != NULL;
    carquet_status_t st = CARQUET_OK;
    if (ok) {
        (void)carquet_writer_write_batch(w, 0, d.ids, ROWS, NULL, NULL);
        (void)carquet_writer_write_batch(w, 1, d.values, ROWS, NULL, NULL);
        (void)carquet_writer_new_row_group(w);
        (void)carquet_writer_write_batch(w, 0, d.ids, 1000, NULL, NULL);
        (void)carquet_writer_write_batch(w, 1, d.values, 1000, NULL, NULL);
        st = carquet_writer_close(w);
    }
    if (s) carquet_schema_free(s);
    fclose(full);
    data_free(&d);
    if (!ok) TEST_FAIL("async_write_error", "setup failed");
    if (st == CARQUET_OK) TEST_FAIL("async_write_error", "close() hid the write failure");
    TEST_PASS("async_write_error");
    return 0;
#endif
}

int main(void) {
    int failures = 0;
    cfg_t base = {CARQUET_COMPRESSION_LZ4_RAW, 0, 0, 1, true, 0};
    failures += test_page_parallel_identical(&base, "page_parallel_lz4");
    cfg_t z = base; z.codec = CARQUET_COMPRESSION_ZSTD;
    failures += test_page_parallel_identical(&z, "page_parallel_zstd");
    cfg_t sn = base; sn.codec = CARQUET_COMPRESSION_SNAPPY; sn.write_batch_size = 50000;
    failures += test_page_parallel_identical(&sn, "page_parallel_snappy_batch50k");
    cfg_t rows = base; rows.max_rows_per_page = 70000;
    failures += test_page_parallel_identical(&rows, "page_parallel_max_rows_70k");
    cfg_t v2 = base; v2.data_page_version = 2; v2.codec = CARQUET_COMPRESSION_GZIP;
    failures += test_page_parallel_identical(&v2, "page_parallel_gzip_v2");
    cfg_t unc = base; unc.codec = CARQUET_COMPRESSION_UNCOMPRESSED;
    failures += test_page_parallel_identical(&unc, "page_parallel_uncompressed");
    failures += test_async_io_identical(CARQUET_COMPRESSION_ZSTD, "async_io_identical_zstd");
    failures += test_async_io_identical(CARQUET_COMPRESSION_UNCOMPRESSED, "async_io_identical_none");
    failures += test_async_buffer_writer();
    failures += test_mixed_nullable();
    failures += test_async_write_error();
    if (failures) { printf("\n%d test(s) FAILED\n", failures); return 1; }
    printf("\nAll parallel-write tests passed\n");
    return 0;
}
