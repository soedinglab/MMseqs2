/**
 * @file fuzz_reader_file.c
 * @brief Fuzz target for the path / FILE* / mmap reader entry points.
 *
 * fuzz_reader opens its input with carquet_reader_open_buffer(), which serves
 * every page from the in-memory (mmap-style) branch of page_reader.c. The
 * default real-world path — carquet_reader_open(path) without mmap — reads
 * page headers and payloads with fseek/fread through a separate set of
 * functions (read_and_parse_page_header_fread() and friends) that fuzz_reader
 * never reaches. This target writes each input to a temp file and drives:
 *
 *   - carquet_get_file_info() / carquet_validate_file()  (path-only helpers)
 *   - carquet_reader_open_file(FILE*) with use_mmap = false: low-level column
 *     reads incl. carquet_column_skip / has_next / remaining, the row-group
 *     statistics API (column_statistics, row_group_matches,
 *     filter_row_groups), geospatial statistics, Arrow per-field metadata and
 *     type refinements, zero-copy probing, bloom-filter float probes,
 *     list_offsets / count_rows over decoded repetition levels, and
 *     carquet_reader_read_arrow() on the fread path
 *   - carquet_reader_open(path) with use_mmap = true: batch reader incl.
 *     carquet_batch_reader_last_error() after a failed next()
 *
 * The whole input is the file, so real Parquet files are valid seeds. Knobs
 * (mmap batch size, thread count, checksum verification) are derived from the
 * input length rather than from file bytes so they never perturb the file.
 */

#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>
#include <carquet/carquet.h>

#define MAX_COLS 32
#define MAX_RGS 4
#define MAX_ROWS 4096

static char g_path[128];

static int write_input(const uint8_t* data, size_t size) {
    if (!g_path[0]) {
        snprintf(g_path, sizeof(g_path), "/tmp/carquet_fuzz_reader_file_%ld.parquet",
                 (long)getpid());
    }
    FILE* f = fopen(g_path, "wb");
    if (!f) return -1;
    size_t w = fwrite(data, 1, size, f);
    fclose(f);
    return w == size ? 0 : -1;
}

/* Predicate value of the right width for the leaf's physical type. */
static int32_t predicate_value(carquet_physical_type_t t, const uint8_t* src,
                               size_t src_len, uint8_t out[16]) {
    memset(out, 0, 16);
    memcpy(out, src, src_len < 16 ? src_len : 16);
    switch (t) {
        case CARQUET_PHYSICAL_BOOLEAN: out[0] &= 1; return 1;
        case CARQUET_PHYSICAL_INT32:
        case CARQUET_PHYSICAL_FLOAT:   return 4;
        case CARQUET_PHYSICAL_INT64:
        case CARQUET_PHYSICAL_DOUBLE:  return 8;
        case CARQUET_PHYSICAL_INT96:   return 12;
        default:                       return (int32_t)(src_len < 16 ? src_len : 16);
    }
}

static void drive_fread_reader(carquet_reader_t* reader, const uint8_t* data,
                               size_t size) {
    const carquet_schema_t* schema = carquet_reader_schema(reader);
    int32_t num_cols = carquet_reader_num_columns(reader);
    int32_t num_rgs = carquet_reader_num_row_groups(reader);
    (void)carquet_reader_is_mmap(reader);

    int16_t* def = malloc(MAX_ROWS * sizeof(int16_t));
    int16_t* rep = malloc(MAX_ROWS * sizeof(int16_t));
    int64_t* offs = malloc((MAX_ROWS + 1) * sizeof(int64_t));
    uint8_t* values = malloc((size_t)MAX_ROWS * 16);
    if (!def || !rep || !offs || !values) goto done;

    for (int32_t c = 0; c < num_cols && c < MAX_COLS; c++) {
        if (schema) {
            (void)carquet_schema_column_name(schema, c);
            (void)carquet_schema_max_def_level(schema, c);
            (void)carquet_schema_max_rep_level(schema, c);
        }
        int32_t nmd = carquet_reader_column_num_metadata(reader, c);
        for (int32_t i = 0; i < nmd && i < 8; i++) {
            const char *k, *v;
            if (carquet_reader_column_get_metadata(reader, c, i, &k, &v) == CARQUET_OK) {
                (void)carquet_reader_column_find_metadata(reader, c, k);
            }
        }
        (void)carquet_reader_column_find_metadata(reader, c, "Label");
        (void)carquet_reader_column_arrow_type_refinement(reader, c);
    }

    for (int32_t rg = 0; rg < num_rgs && rg < MAX_RGS; rg++) {
        for (int32_t c = 0; c < num_cols && c < MAX_COLS; c++) {
            (void)carquet_reader_can_zero_copy(reader, rg, c);

            carquet_column_statistics_t st;
            if (carquet_reader_column_statistics(reader, rg, c, &st) == CARQUET_OK &&
                st.has_min_max && st.min_value && st.min_value_size > 0) {
                /* Volatile read of the stats bytes: they point into
                 * reader-owned metadata and must be fully readable. */
                volatile uint8_t sink = 0;
                for (int32_t i = 0; i < st.min_value_size; i++)
                    sink ^= ((const uint8_t*)st.min_value)[i];
                for (int32_t i = 0; st.max_value && i < st.max_value_size; i++)
                    sink ^= ((const uint8_t*)st.max_value)[i];
                (void)sink;
            }

            carquet_physical_type_t t = schema ? carquet_schema_column_type(schema, c)
                                               : CARQUET_PHYSICAL_INT32;
            uint8_t pv[16];
            int32_t pv_len = predicate_value(t, data + (size_t)c % (size - 8), 8, pv);
            bool might = false;
            (void)carquet_reader_row_group_matches(
                reader, rg, c, (carquet_compare_op_t)((size + (size_t)c) % 6),
                pv, pv_len, &might);
            if (rg == 0) {
                int32_t idx[MAX_RGS * 2];
                int32_t n = carquet_reader_filter_row_groups(
                    reader, c, CARQUET_COMPARE_GE, pv, pv_len, idx, MAX_RGS * 2);
                assert(n <= MAX_RGS * 2);
                for (int32_t i = 0; i < n; i++) assert(idx[i] >= 0 && idx[i] < num_rgs);
            }

            carquet_geospatial_statistics_t geo;
            (void)carquet_reader_geospatial_statistics(reader, rg, c, &geo);

            carquet_bloom_filter_t* bf = carquet_reader_get_bloom_filter(reader, rg, c, NULL);
            if (bf) {
                (void)carquet_bloom_filter_check_float(bf, 1.5f);
                (void)carquet_bloom_filter_check_i32(bf, 7);
                carquet_bloom_filter_destroy(bf);
            }

            carquet_error_t err = CARQUET_ERROR_INIT;
            carquet_column_reader_t* cr = carquet_reader_get_column(reader, rg, c, &err);
            if (!cr) continue;

            /* Skip a prefix, then read the rest in two batches. */
            (void)carquet_column_remaining(cr);
            int64_t skip_n = (int64_t)(size % 17);
            int64_t skipped = carquet_column_skip(cr, skip_n);
            assert(skipped >= 0 && skipped <= skip_n);
            int64_t total = 0;
            for (int pass = 0; pass < 2 && carquet_column_has_next(cr); pass++) {
                carquet_error_t rerr = CARQUET_ERROR_INIT;
                int64_t n = carquet_column_read_batch_ex(cr, values, MAX_ROWS / 2,
                                                         def, rep, &rerr);
                if (n < 0) { assert(rerr.code != CARQUET_OK); break; }
                if (n == 0) break;
                assert(n <= MAX_ROWS / 2);
                if (schema && carquet_schema_max_rep_level(schema, c) > 0) {
                    (void)carquet_count_rows(rep, n);
                    int64_t nl = carquet_list_offsets(rep, n, 1, offs, MAX_ROWS + 1);
                    assert(nl >= 0 && nl <= MAX_ROWS + 1);
                }
                total += n;
            }
            (void)total;
            carquet_column_reader_free(cr);
        }
    }

    /* Nested Arrow reassembly on the fread path (first row group only). */
    if (num_rgs > 0) {
        struct ArrowSchema as;
        struct ArrowArray aa;
        memset(&as, 0, sizeof(as));
        memset(&aa, 0, sizeof(aa));
        if (carquet_reader_read_arrow(reader, 0, &as, &aa, NULL) == CARQUET_OK) {
            if (aa.release) aa.release(&aa);
            if (as.release) as.release(&as);
        }
    }

done:
    free(def);
    free(rep);
    free(offs);
    free(values);
}

static void drive_mmap_reader(carquet_reader_t* reader, size_t size) {
    int32_t num_cols = carquet_reader_num_columns(reader);
    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    cfg.batch_size = (int32_t)(1 + (size % 7) * 97);
    cfg.num_threads = (int32_t)(size % 3);
    cfg.use_mmap = true;

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_batch_reader_t* br = carquet_batch_reader_create(reader, &cfg, &err);
    if (!br) return;
    carquet_row_batch_t* batch = NULL;
    carquet_status_t st = CARQUET_OK;
    for (int i = 0; i < 8; i++) {
        st = carquet_batch_reader_next(br, &batch);
        if (st != CARQUET_OK || !batch) break;
        for (int32_t c = 0; c < num_cols && c < MAX_COLS; c++) {
            const void* d; const uint8_t* nulls; int64_t count;
            if (carquet_row_batch_column(batch, c, &d, &nulls, &count) == CARQUET_OK) {
                assert(count >= 0);
            }
        }
        (void)carquet_row_batch_num_columns(batch);
        carquet_row_batch_free(batch);
        batch = NULL;
    }
    if (st != CARQUET_OK && st != CARQUET_ERROR_END_OF_DATA) {
        const carquet_error_t* le = carquet_batch_reader_last_error(br);
        if (le) {
            char buf[256];
            (void)carquet_error_format(le, buf, sizeof(buf));
        }
    }
    carquet_batch_reader_free(br);
}

int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
    if (size < 12) return 0;
    (void)carquet_init();
    if (write_input(data, size) != 0) return 0;

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_file_info_t info;
    (void)carquet_get_file_info(g_path, &info, &err);
    err = (carquet_error_t)CARQUET_ERROR_INIT;
    (void)carquet_validate_file(g_path, &err);

    carquet_reader_options_t opts;
    carquet_reader_options_init(&opts);
    opts.use_mmap = false;
    opts.verify_checksums = (size & 1) != 0;
    opts.num_threads = 1;

    FILE* f = fopen(g_path, "rb");
    if (f) {
        err = (carquet_error_t)CARQUET_ERROR_INIT;
        carquet_reader_t* r = carquet_reader_open_file(f, &opts, &err);
        if (r) {
            drive_fread_reader(r, data, size);
            carquet_reader_close(r);
        } else {
            char buf[256];
            (void)carquet_error_format(&err, buf, sizeof(buf));
        }
        fclose(f);
    }

    opts.use_mmap = true;
    err = (carquet_error_t)CARQUET_ERROR_INIT;
    carquet_reader_t* r = carquet_reader_open(g_path, &opts, &err);
    if (r) {
        drive_mmap_reader(r, size);
        carquet_reader_close(r);
    }
    return 0;
}

#ifdef AFL_MAIN
#include <sys/stat.h>
int main(int argc, char** argv) {
    if (argc != 2) { fprintf(stderr, "Usage: %s <input_file>\n", argv[0]); return 1; }
    FILE* f = fopen(argv[1], "rb");
    if (!f) { perror("fopen"); return 1; }
    struct stat st; fstat(fileno(f), &st);
    uint8_t* d = malloc((size_t)st.st_size);
    if (!d) { fclose(f); return 1; }
    fread(d, 1, (size_t)st.st_size, f); fclose(f);
    int r = LLVMFuzzerTestOneInput(d, (size_t)st.st_size);
    free(d); return r;
}
#endif
