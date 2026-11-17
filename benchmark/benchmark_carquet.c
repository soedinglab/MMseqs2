/**
 * @file benchmark_carquet.c
 * @brief Performance benchmarks for Carquet
 *
 * Methodology:
 *   - Write and read are benchmarked separately
 *   - Page cache is purged between write and read phases (macOS: F_NOCACHE,
 *     Windows: unbuffered copy, see bench_purge_win32.h)
 *   - 3 warmup + 11 measured iterations, trimmed mean (drop min + max, median of 9)
 *   - Data is pre-generated outside timing loops
 *   - Reads are timed twice: with use_mmap = true (the "read" column, as
 *     before) and with use_mmap = false (the C default, "read_nommap")
 *
 * Nullable datasets ("<size>_nullable"): same columns and values, but OPTIONAL
 * with ~NULL_PERCENT % nulls per column (deterministic hash of the row index,
 * shared with the C++/Python benchmarks) and RLE_DICTIONARY on `category`.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#ifdef _WIN32
#include <windows.h>
#include <sys/stat.h>
#include "bench_purge_win32.h"
#else
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#endif

#ifdef __APPLE__
#include <sys/sysctl.h>
#endif

#include <carquet/carquet.h>

#define WARMUP_ITERATIONS 3
#define BENCH_ITERATIONS_SMALL 51   /* sub-ms ops need more samples */
#define BENCH_ITERATIONS_MEDIUM 21
#define BENCH_ITERATIONS_LARGE 11
#define MAX_BENCH_ITERATIONS 201

static int compare_double(const void* a, const void* b) {
    double da = *(const double*)a;
    double db = *(const double*)b;
    return (da > db) - (da < db);
}

/* Trimmed median: sort, drop min and max, return median of the rest */
static double trimmed_median(double* arr, int n) {
    qsort(arr, n, sizeof(double), compare_double);
    if (n <= 2) return arr[n / 2];
    /* Trim 1 from each end */
    int lo = 1, hi = n - 2;
    int trimmed = hi - lo + 1;
    if (trimmed % 2 == 1) return arr[lo + trimmed / 2];
    return (arr[lo + trimmed / 2 - 1] + arr[lo + trimmed / 2]) / 2.0;
}

typedef struct {
    const char* name;
    int rows;
    int nullable;
} dataset_t;

typedef struct {
    carquet_compression_t codec;
    const char* name;
} compression_config_t;

static int get_benchmark_zstd_level(void) {
    const char* env = getenv("CARQUET_BENCH_ZSTD_LEVEL");
    if (!env || env[0] == '\0') {
        return 1;
    }

    char* end = NULL;
    long level = strtol(env, &end, 10);
    if (end == env || *end != '\0') {
        return 1;
    }
    if (level < 1) level = 1;
    if (level > 22) level = 22;
    return (int)level;
}

static double get_time_ms(void) {
#ifdef _WIN32
    static LARGE_INTEGER freq = {0};
    LARGE_INTEGER counter;
    if (freq.QuadPart == 0) {
        QueryPerformanceFrequency(&freq);
    }
    QueryPerformanceCounter(&counter);
    return (double)counter.QuadPart * 1000.0 / (double)freq.QuadPart;
#else
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec * 1000.0 + ts.tv_nsec / 1000000.0;
#endif
}

static long get_file_size(const char* filename) {
    struct stat st;
    if (stat(filename, &st) == 0) {
        return (long)st.st_size;
    }
    return 0;
}

static const char* get_temp_dir(void) {
    const char* override = getenv("CARQUET_BENCH_TMPDIR");
    if (override && override[0] != '\0') {
        return override;
    }
#ifdef _WIN32
    static char temp_dir[512] = {0};
    if (temp_dir[0] == 0) {
        const char* tmp = getenv("TEMP");
        if (!tmp) tmp = getenv("TMP");
        if (!tmp) tmp = ".";
        snprintf(temp_dir, sizeof(temp_dir), "%s", tmp);
    }
    return temp_dir;
#else
    return "/tmp";
#endif
}

/**
 * Evict a file from the OS page cache so reads hit storage, not DRAM.
 *
 * macOS: No unprivileged API reliably evicts unified buffer cache pages.
 *        We copy the file to a new path with F_NOCACHE, delete the original,
 *        and rename. The new fd's pages won't be in cache.
 * Linux: posix_fadvise DONTNEED.
 * Windows: same copy-and-replace as macOS with unbuffered I/O (bench_purge_win32.h).
 */
static void purge_file_cache(const char* filename) {
#if defined(__APPLE__)
    /* Strategy: copy file bypassing cache, replace original */
    char tmp[512];
    snprintf(tmp, sizeof(tmp), "%s.nocache", filename);

    int src = open(filename, O_RDONLY);
    if (src < 0) return;

    struct stat st;
    fstat(src, &st);

    int dst = open(tmp, O_WRONLY | O_CREAT | O_TRUNC, 0644);
    if (dst < 0) { close(src); return; }
    fcntl(dst, F_NOCACHE, 1);  /* Don't cache the written pages */
    fcntl(src, F_NOCACHE, 1);  /* Hint to evict source pages as we read */

    char buf[262144];
    ssize_t n;
    while ((n = read(src, buf, sizeof(buf))) > 0) {
        write(dst, buf, n);
    }
    close(src);
    close(dst);

    /* Replace original with the uncached copy */
    unlink(filename);
    rename(tmp, filename);
#elif defined(__linux__)
    int fd = open(filename, O_RDONLY);
    if (fd < 0) return;
    struct stat st;
    if (fstat(fd, &st) == 0) {
        posix_fadvise(fd, 0, st.st_size, POSIX_FADV_DONTNEED);
    }
    close(fd);
#elif defined(_WIN32)
    bench_purge_file_cache_win32(filename);
#else
    (void)filename;
#endif
}

/* Simple LCG random for reproducible results */
static uint32_t lcg_state = 42;
static void lcg_seed(uint32_t seed) { lcg_state = seed; }
static uint32_t lcg_rand(void) {
    lcg_state = lcg_state * 1103515245 + 12345;
    return (lcg_state >> 16) & 0x7FFF;
}
static double lcg_normal(void) {
    double u1 = (lcg_rand() + 1.0) / 32768.0;
    double u2 = (lcg_rand() + 1.0) / 32768.0;
    return sqrt(-2.0 * log(u1)) * cos(2.0 * 3.14159265358979 * u2);
}

/* Null pattern of the nullable datasets. Must match bench_common.null_masks()
 * and benchmark_arrow_cpp.cpp: lowbias32 hash of (row ^ per-column seed). */
#define NULL_PERCENT 7
#define NUM_COLUMNS 3

static int bench_is_null(uint32_t row, int col) {
    uint32_t x = row ^ (0x9E3779B9u * (uint32_t)(col + 1));
    x ^= x >> 16;
    x *= 0x7FEB352Du;
    x ^= x >> 15;
    x *= 0x846CA68Bu;
    x ^= x >> 16;
    return (x % 100u) < NULL_PERCENT;
}

static int64_t bench_rg_rows(int num_rows) {
    return num_rows > 1000000 ? num_rows / 10 : 100000;
}

/* Pre-generated test data (allocated once, reused across iterations) */
typedef struct {
    int64_t* ids;
    double* values;
    int32_t* categories;
    int num_rows;
    /* Nullable datasets only: value arrays above are compacted to the
     * non-null values (carquet's write_batch input), def_levels hold one
     * entry per row and rg_value_start[c][rg] is the first value of row
     * group rg in column c. All built outside the timed loops. */
    int nullable;
    int16_t* def_levels[NUM_COLUMNS];
    int64_t* rg_value_start[NUM_COLUMNS];
    int64_t null_counts[NUM_COLUMNS];
} test_data_t;

static void compact_nullable_column(test_data_t* td, int col, void* data, size_t elem_size) {
    int64_t rg_rows = bench_rg_rows(td->num_rows);
    int64_t num_rgs = (td->num_rows + rg_rows - 1) / rg_rows;
    int16_t* def = malloc((size_t)td->num_rows * sizeof(int16_t));
    int64_t* rg_start = malloc((size_t)(num_rgs + 1) * sizeof(int64_t));
    uint8_t* bytes = (uint8_t*)data;
    int64_t out = 0;

    for (int64_t i = 0; i < td->num_rows; i++) {
        if (i % rg_rows == 0) {
            rg_start[i / rg_rows] = out;
        }
        if (bench_is_null((uint32_t)i, col)) {
            def[i] = 0;
        } else {
            def[i] = 1;
            if (out != i) {
                memcpy(bytes + out * elem_size, bytes + i * elem_size, elem_size);
            }
            out++;
        }
    }
    rg_start[num_rgs] = out;

    td->def_levels[col] = def;
    td->rg_value_start[col] = rg_start;
    td->null_counts[col] = td->num_rows - out;
}

static test_data_t* test_data_create(int num_rows, int nullable) {
    test_data_t* td = calloc(1, sizeof(test_data_t));
    td->ids = malloc(num_rows * sizeof(int64_t));
    td->values = malloc(num_rows * sizeof(double));
    td->categories = malloc(num_rows * sizeof(int32_t));
    td->num_rows = num_rows;
    td->nullable = nullable;

    lcg_seed(42);
    for (int i = 0; i < num_rows; i++) {
        td->ids[i] = 1000000 + (lcg_rand() % 9000000);
        td->values[i] = fabs(100.0 + 50.0 * lcg_normal());
        td->categories[i] = lcg_rand() % 100;
    }

    if (nullable) {
        compact_nullable_column(td, 0, td->ids, sizeof(int64_t));
        compact_nullable_column(td, 1, td->values, sizeof(double));
        compact_nullable_column(td, 2, td->categories, sizeof(int32_t));
    }
    return td;
}

static void test_data_destroy(test_data_t* td) {
    if (td) {
        free(td->ids);
        free(td->values);
        free(td->categories);
        for (int c = 0; c < NUM_COLUMNS; c++) {
            free(td->def_levels[c]);
            free(td->rg_value_start[c]);
        }
        free(td);
    }
}

/* Portable popcount (no __builtin_popcountll on MSVC). */
static int popcount64(uint64_t x) {
    x = x - ((x >> 1) & 0x5555555555555555ULL);
    x = (x & 0x3333333333333333ULL) + ((x >> 2) & 0x3333333333333333ULL);
    x = (x + (x >> 4)) & 0x0F0F0F0F0F0F0F0FULL;
    return (int)((x * 0x0101010101010101ULL) >> 56);
}

/* Nulls in the first `count` bits of a validity bitmap (bit set = present). */
static int64_t count_nulls(const uint8_t* bitmap, int64_t count) {
    if (!bitmap || count <= 0) return 0;
    int64_t valid = 0;
    int64_t full_bytes = count / 8;
    int64_t i = 0;
    for (; i + 8 <= full_bytes; i += 8) {
        uint64_t word;
        memcpy(&word, bitmap + i, sizeof(word));
        valid += popcount64(word);
    }
    for (; i < full_bytes; i++) {
        valid += popcount64(bitmap[i]);
    }
    int rem = (int)(count % 8);
    if (rem) {
        valid += popcount64(bitmap[full_bytes] & ((1u << rem) - 1u));
    }
    return count - valid;
}

static double benchmark_write(const char* filename, const test_data_t* td,
                               carquet_compression_t codec) {
    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_field_repetition_t rep = td->nullable ? CARQUET_REPETITION_OPTIONAL
                                                  : CARQUET_REPETITION_REQUIRED;
    carquet_schema_t* schema = carquet_schema_create(&err);
    (void)carquet_schema_add_column(schema, "id", CARQUET_PHYSICAL_INT64,
                                    NULL, rep, 0, 0);
    (void)carquet_schema_add_column(schema, "value", CARQUET_PHYSICAL_DOUBLE,
                                    NULL, rep, 0, 0);
    (void)carquet_schema_add_column(schema, "category", CARQUET_PHYSICAL_INT32,
                                    NULL, rep, 0, 0);

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = codec;
    int64_t rg_rows = td->num_rows > 1000000 ? td->num_rows / 10 : 100000;
    if (codec == CARQUET_COMPRESSION_ZSTD) {
        opts.compression_level = get_benchmark_zstd_level();
    }
    opts.row_group_size = rg_rows * (int64_t)(sizeof(int64_t) + sizeof(double) + sizeof(int32_t));

    double start = get_time_ms();

    carquet_writer_t* writer = carquet_writer_create(filename, schema, &opts, &err);
    if (!writer) {
        fprintf(stderr, "carquet_writer_create failed for %s: %s\n",
                filename, err.message[0] ? err.message : "unknown error");
        carquet_schema_free(schema);
        return 0.0;
    }
    if (td->nullable &&
        carquet_writer_set_column_encoding(writer, 2, CARQUET_ENCODING_RLE_DICTIONARY) != CARQUET_OK) {
        fprintf(stderr, "carquet_writer_set_column_encoding failed for %s\n", filename);
        (void)carquet_writer_close(writer);
        carquet_schema_free(schema);
        return 0.0;
    }
    for (int64_t offset = 0; offset < td->num_rows; offset += rg_rows) {
        int64_t chunk = td->num_rows - offset;
        if (chunk > rg_rows) {
            chunk = rg_rows;
        }

        int ok;
        if (td->nullable) {
            int64_t rg = offset / rg_rows;
            ok = carquet_writer_write_batch(writer, 0, td->ids + td->rg_value_start[0][rg], chunk,
                                            td->def_levels[0] + offset, NULL) == CARQUET_OK &&
                 carquet_writer_write_batch(writer, 1, td->values + td->rg_value_start[1][rg], chunk,
                                            td->def_levels[1] + offset, NULL) == CARQUET_OK &&
                 carquet_writer_write_batch(writer, 2, td->categories + td->rg_value_start[2][rg], chunk,
                                            td->def_levels[2] + offset, NULL) == CARQUET_OK;
        } else {
            ok = carquet_writer_write_batch(writer, 0, td->ids + offset, chunk, NULL, NULL) == CARQUET_OK &&
                 carquet_writer_write_batch(writer, 1, td->values + offset, chunk, NULL, NULL) == CARQUET_OK &&
                 carquet_writer_write_batch(writer, 2, td->categories + offset, chunk, NULL, NULL) == CARQUET_OK;
        }
        if (!ok) {
            fprintf(stderr, "carquet_writer_write_batch failed for %s: %s\n",
                    filename, err.message[0] ? err.message : "unknown error");
            (void)carquet_writer_close(writer);
            carquet_schema_free(schema);
            return 0.0;
        }

        if (offset + chunk < td->num_rows) {
            if (carquet_writer_new_row_group(writer) != CARQUET_OK) {
                fprintf(stderr, "carquet_writer_new_row_group failed for %s: %s\n",
                        filename, err.message[0] ? err.message : "unknown error");
                (void)carquet_writer_close(writer);
                carquet_schema_free(schema);
                return 0.0;
            }
        }
    }
    if (carquet_writer_close(writer) != CARQUET_OK) {
        fprintf(stderr, "carquet_writer_close failed for %s: %s\n",
                filename, err.message[0] ? err.message : "unknown error");
        carquet_schema_free(schema);
        return 0.0;
    }

    double elapsed = get_time_ms() - start;
    carquet_schema_free(schema);
    return elapsed;
}

/* Returns the elapsed ms, 0.0 if the file could not be opened, or -1.0 if
 * expected_nulls is given and the per-column null counts do not match. */
static double benchmark_read(const char* filename, int expected_rows,
                              carquet_thread_pool_t* pool, bool use_mmap,
                              const int64_t* expected_nulls) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    int null_mismatch = 0;

    double start = get_time_ms();

    carquet_reader_options_t opts;
    carquet_reader_options_init(&opts);
    opts.use_mmap = use_mmap;
    opts.verify_checksums = false;

    carquet_reader_t* reader = carquet_reader_open(filename, &opts, &err);
    if (!reader) {
        fprintf(stderr, "carquet_reader_open failed for %s: %s\n",
                filename, err.message[0] ? err.message : "unknown error");
        return 0.0;
    }

    carquet_batch_reader_config_t config;
    carquet_batch_reader_config_init(&config);
    config.batch_size = 262144;
    config.thread_pool = pool;

    carquet_batch_reader_t* batch_reader = carquet_batch_reader_create(reader, &config, &err);
    if (batch_reader) {
        carquet_row_batch_t* batch = NULL;
        int64_t total_rows = 0;
        int64_t nulls_seen[NUM_COLUMNS] = {0, 0, 0};
        volatile int64_t checksum = 0;

        while (carquet_batch_reader_next(batch_reader, &batch) == CARQUET_OK && batch) {
            int64_t batch_rows = carquet_row_batch_num_rows(batch);
            total_rows += batch_rows;

            const void* data;
            const uint8_t* nulls;
            int64_t count;

            /* nulls is NULL for REQUIRED columns, so count_nulls() only
             * does work (one popcount pass over the bitmap) on OPTIONAL ones. */
            if (carquet_row_batch_column(batch, 0, &data, &nulls, &count) == CARQUET_OK && data) {
                const int64_t* ids = (const int64_t*)data;
                checksum += ids[0] + ids[count > 1 ? count - 1 : 0];
                nulls_seen[0] += count_nulls(nulls, count);
            }
            if (carquet_row_batch_column(batch, 1, &data, &nulls, &count) == CARQUET_OK && data) {
                const double* values = (const double*)data;
                checksum += (int64_t)values[0] + (int64_t)values[count > 1 ? count - 1 : 0];
                nulls_seen[1] += count_nulls(nulls, count);
            }
            if (carquet_row_batch_column(batch, 2, &data, &nulls, &count) == CARQUET_OK && data) {
                const int32_t* cats = (const int32_t*)data;
                checksum += cats[0] + cats[count > 1 ? count - 1 : 0];
                nulls_seen[2] += count_nulls(nulls, count);
            }

            carquet_row_batch_free(batch);
            batch = NULL;
        }

        checksum += nulls_seen[0] + nulls_seen[1] + nulls_seen[2];
        (void)checksum;
        if (total_rows != expected_rows) {
            fprintf(stderr, "Warning: row count mismatch %lld vs %d\n",
                    (long long)total_rows, expected_rows);
        }
        if (expected_nulls) {
            for (int c = 0; c < NUM_COLUMNS; c++) {
                if (nulls_seen[c] != expected_nulls[c]) {
                    fprintf(stderr, "Null count mismatch in column %d of %s: %lld vs %lld\n",
                            c, filename, (long long)nulls_seen[c], (long long)expected_nulls[c]);
                    null_mismatch = 1;
                }
            }
        }
        carquet_batch_reader_free(batch_reader);
    }

    carquet_reader_close(reader);
    double elapsed = get_time_ms() - start;
    return null_mismatch ? -1.0 : elapsed;
}

static void run_benchmark(const char* dataset_name, int num_rows, int nullable,
                          carquet_compression_t codec, const char* compression_name) {
    char filename[256];
    snprintf(filename, sizeof(filename), "%s/benchmark_%s_%s_carquet.parquet",
             get_temp_dir(), dataset_name, compression_name);

    printf("\n=== %s (%d rows, %s) ===\n", dataset_name, num_rows, compression_name);

    test_data_t* td = test_data_create(num_rows, nullable);
    const int64_t* expected_nulls = nullable ? td->null_counts : NULL;

    /* Scale iterations: small datasets need more samples for stability.
     * CARQUET_BENCH_ITERATIONS env var overrides all sizes. */
    int iters;
    const char* iter_env = getenv("CARQUET_BENCH_ITERATIONS");
    if (iter_env && atoi(iter_env) > 0) {
        iters = atoi(iter_env);
        if (iters > MAX_BENCH_ITERATIONS) iters = MAX_BENCH_ITERATIONS;
    } else if (num_rows <= 100000) {
        iters = BENCH_ITERATIONS_SMALL;
    } else if (num_rows <= 1000000) {
        iters = BENCH_ITERATIONS_MEDIUM;
    } else {
        iters = BENCH_ITERATIONS_LARGE;  /* large and xlarge */
    }

    /* Warmup writes */
    for (int i = 0; i < WARMUP_ITERATIONS; i++) {
        benchmark_write(filename, td, codec);
    }

    /* Benchmark writes */
    double write_times[MAX_BENCH_ITERATIONS];
    for (int i = 0; i < iters; i++) {
        write_times[i] = benchmark_write(filename, td, codec);
    }

    long file_size = get_file_size(filename);

    /* Create a reusable thread pool — avoids pthread create/join per iteration */
    carquet_thread_pool_t* pool = carquet_thread_pool_create(0);

    /* Purge cache once, then warmup reads (first read = cold, rest = warm) */
    purge_file_cache(filename);

    int read_failed = 0;
    for (int i = 0; i < WARMUP_ITERATIONS; i++) {
        if (benchmark_read(filename, num_rows, pool, true, expected_nulls) < 0.0) read_failed = 1;
    }

    /* Benchmark reads (warm cache — realistic for most workloads) */
    double read_times[MAX_BENCH_ITERATIONS];
    for (int i = 0; i < iters; i++) {
        read_times[i] = benchmark_read(filename, num_rows, pool, true, expected_nulls);
        if (read_times[i] < 0.0) read_failed = 1;
    }

    /* Same protocol without mmap (use_mmap = false is the C default): purge,
     * warm up, measure. Runs after the mmap pass so that one is unchanged. */
    purge_file_cache(filename);

    for (int i = 0; i < WARMUP_ITERATIONS; i++) {
        if (benchmark_read(filename, num_rows, pool, false, expected_nulls) < 0.0) read_failed = 1;
    }

    double read_nommap_times[MAX_BENCH_ITERATIONS];
    for (int i = 0; i < iters; i++) {
        read_nommap_times[i] = benchmark_read(filename, num_rows, pool, false, expected_nulls);
        if (read_nommap_times[i] < 0.0) read_failed = 1;
    }

    carquet_thread_pool_destroy(pool);

    if (read_failed) {
        fprintf(stderr, "Read validation failed for %s\n", filename);
        test_data_destroy(td);
        return;
    }

    double write_med = trimmed_median(write_times, iters);
    double read_med = trimmed_median(read_times, iters);
    double read_nommap_med = trimmed_median(read_nommap_times, iters);

    double rows_per_sec_write = (num_rows / write_med) * 1000;
    double rows_per_sec_read = (num_rows / read_med) * 1000;

    printf("  Write: %.2f ms (%.2f M rows/sec)\n", write_med, rows_per_sec_write / 1e6);
    printf("  Read:  %.2f ms (%.2f M rows/sec)\n", read_med, rows_per_sec_read / 1e6);
    printf("  Read (no mmap): %.2f ms (%.2f M rows/sec)\n",
           read_nommap_med, (num_rows / read_nommap_med) * 1000 / 1e6);
    printf("  File:  %.2f MB (%.2f bytes/row)\n",
           file_size / (1024.0 * 1024.0), (double)file_size / num_rows);

    /* 8th field (read_nommap) is optional for the orchestrator's parser. */
    printf("CSV:carquet,%s,%s,%d,%.2f,%.2f,%ld,%.2f\n",
           dataset_name, compression_name, num_rows, write_med, read_med, file_size,
           read_nommap_med);

    /* File is NOT removed here — the orchestrator handles cleanup
     * after cross-read benchmarking is done. */
    test_data_destroy(td);
}

static int find_dataset(const char* name, dataset_t* datasets, int n) {
    for (int i = 0; i < n; i++)
        if (strcmp(datasets[i].name, name) == 0) return i;
    return -1;
}

static int find_compression(const char* name, compression_config_t* comps, int n) {
    for (int i = 0; i < n; i++)
        if (strcmp(comps[i].name, name) == 0) return i;
    return -1;
}

static int run_crossread(const char* filename) {
    carquet_error_t err = CARQUET_ERROR_INIT;

    /* Open file to get row count from metadata */
    carquet_reader_options_t opts;
    carquet_reader_options_init(&opts);
    opts.use_mmap = true;
    opts.verify_checksums = false;

    carquet_reader_t* reader = carquet_reader_open(filename, &opts, &err);
    if (!reader) {
        fprintf(stderr, "carquet_reader_open failed for %s: %s\n",
                filename, err.message[0] ? err.message : "unknown error");
        return 1;
    }

    int64_t num_rows = carquet_reader_num_rows(reader);
    carquet_reader_close(reader);

    if (num_rows <= 0) {
        fprintf(stderr, "File %s has no rows\n", filename);
        return 1;
    }

    int expected_rows = (int)num_rows;
    long file_size = get_file_size(filename);

    /* Scale iterations based on row count */
    int iters;
    const char* iter_env = getenv("CARQUET_BENCH_ITERATIONS");
    if (iter_env && atoi(iter_env) > 0) {
        iters = atoi(iter_env);
        if (iters > MAX_BENCH_ITERATIONS) iters = MAX_BENCH_ITERATIONS;
    } else if (num_rows <= 100000) {
        iters = BENCH_ITERATIONS_SMALL;
    } else if (num_rows <= 1000000) {
        iters = BENCH_ITERATIONS_MEDIUM;
    } else {
        iters = BENCH_ITERATIONS_LARGE;
    }

    carquet_thread_pool_t* pool = carquet_thread_pool_create(0);

    /* Purge file cache before reads */
    purge_file_cache(filename);

    /* Warmup reads */
    for (int i = 0; i < WARMUP_ITERATIONS; i++) {
        double t = benchmark_read(filename, expected_rows, pool, true, NULL);
        if (t <= 0.0) {
            fprintf(stderr, "Warmup read failed for %s\n", filename);
            carquet_thread_pool_destroy(pool);
            return 1;
        }
    }

    /* Measured reads */
    double read_times[MAX_BENCH_ITERATIONS];
    for (int i = 0; i < iters; i++) {
        read_times[i] = benchmark_read(filename, expected_rows, pool, true, NULL);
        if (read_times[i] <= 0.0) {
            fprintf(stderr, "Benchmark read failed for %s\n", filename);
            carquet_thread_pool_destroy(pool);
            return 1;
        }
    }

    carquet_thread_pool_destroy(pool);

    double read_med = trimmed_median(read_times, iters);

    printf("CSV:carquet,crossread,external,%d,0.00,%.2f,%ld\n",
           expected_rows, read_med, file_size);

    return 0;
}

int main(int argc, char* argv[]) {
    setvbuf(stdout, NULL, _IONBF, 0);

    if (argc == 2 && strcmp(argv[1], "--version") == 0) {
        printf("%s\n", carquet_version());
        return 0;
    }

    /* read <filename> subcommand for cross-read benchmarking */
    if (argc == 2 && strcmp(argv[1], "read") == 0) {
        fprintf(stderr, "Usage: %s read <filename>\n", argv[0]);
        return 1;
    }
    if (argc >= 3 && strcmp(argv[1], "read") == 0) {
        return run_crossread(argv[2]);
    }

    dataset_t datasets[] = {
        {"small", 100000, 0},
        {"medium", 1000000, 0},
        {"large", 10000000, 0},
        {"xlarge", 100000000, 0},
        /* OPTIONAL columns, ~NULL_PERCENT % nulls, dictionary on category */
        {"small_nullable", 100000, 1},
        {"medium_nullable", 1000000, 1},
        {"large_nullable", 10000000, 1},
        {"xlarge_nullable", 100000000, 1}
    };
    int num_datasets = sizeof(datasets) / sizeof(datasets[0]);

    compression_config_t compressions[] = {
        {CARQUET_COMPRESSION_UNCOMPRESSED, "none"},
        {CARQUET_COMPRESSION_SNAPPY, "snappy"},
        {CARQUET_COMPRESSION_ZSTD, "zstd"},
        {CARQUET_COMPRESSION_LZ4_RAW, "lz4"}
    };
    int num_compressions = sizeof(compressions) / sizeof(compressions[0]);

    if (argc == 3) {
        int d = find_dataset(argv[1], datasets, num_datasets);
        int c = find_compression(argv[2], compressions, num_compressions);
        if (d < 0 || c < 0) {
            fprintf(stderr, "Usage: %s [dataset] [compression]\n", argv[0]);
            fprintf(stderr, "  datasets:     small, medium, large, xlarge (+ _nullable)\n");
            fprintf(stderr, "  compressions: none, snappy, zstd, lz4\n");
            return 1;
        }
        run_benchmark(datasets[d].name, datasets[d].rows, datasets[d].nullable,
                     compressions[c].codec, compressions[c].name);
        return 0;
    }

    printf("Carquet Benchmark\n");

    for (int d = 0; d < num_datasets; d++) {
        if (datasets[d].nullable) continue;  /* opt-in: pass "<size>_nullable <codec>" */
        for (int c = 0; c < num_compressions; c++) {
            run_benchmark(datasets[d].name, datasets[d].rows, datasets[d].nullable,
                         compressions[c].codec, compressions[c].name);
        }
    }

    printf("\nBenchmark complete.\n");
    return 0;
}
