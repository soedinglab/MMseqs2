#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <future>
#include <memory>
#include <string>
#include <thread>
#include <vector>

#ifdef _WIN32
/* windows.h's min/max macros break Arrow's std::numeric_limits<>::max() calls. */
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#include <sys/stat.h>
#include "bench_purge_win32.h"
#else
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

#include <arrow/api.h>
#include <arrow/io/file.h>
#include <arrow/result.h>
#include <arrow/util/config.h>

#include <parquet/arrow/reader.h>
#include <parquet/arrow/writer.h>
#include <parquet/column_reader.h>
#include <parquet/file_reader.h>
#include <parquet/properties.h>

#define WARMUP_ITERATIONS 3
#define BENCH_ITERATIONS_SMALL 51
#define BENCH_ITERATIONS_MEDIUM 21
#define BENCH_ITERATIONS_LARGE 11
#define MAX_BENCH_ITERATIONS 201

/* Null pattern of the nullable datasets. Must match benchmark_carquet.c and
 * bench_common.null_masks(): lowbias32 hash of (row ^ per-column seed). */
#define NULL_PERCENT 7
#define NUM_COLUMNS 3

typedef struct {
    const char* name;
    int rows;
    bool nullable;
} dataset_t;

typedef struct {
    parquet::Compression::type codec;
    const char* name;
} compression_config_t;

typedef struct {
    int num_rows;
    bool nullable;
    int64_t null_counts[NUM_COLUMNS];
    std::shared_ptr<arrow::Table> table;
} test_data_t;

static bool bench_is_null(uint32_t row, int col) {
    uint32_t x = row ^ (0x9E3779B9u * static_cast<uint32_t>(col + 1));
    x ^= x >> 16;
    x *= 0x7FEB352Du;
    x ^= x >> 15;
    x *= 0x846CA68Bu;
    x ^= x >> 16;
    return (x % 100u) < NULL_PERCENT;
}

static int get_benchmark_zstd_level(void) {
    const char* env = std::getenv("CARQUET_BENCH_ZSTD_LEVEL");
    if (!env || env[0] == '\0') {
        return 1;
    }

    char* end = NULL;
    long level = std::strtol(env, &end, 10);
    if (end == env || *end != '\0') {
        return 1;
    }
    if (level < 1) level = 1;
    if (level > 22) level = 22;
    return static_cast<int>(level);
}

static double get_time_ms(void) {
#ifdef _WIN32
    static LARGE_INTEGER freq = {0};
    LARGE_INTEGER counter;
    if (freq.QuadPart == 0) {
        QueryPerformanceFrequency(&freq);
    }
    QueryPerformanceCounter(&counter);
    return static_cast<double>(counter.QuadPart) * 1000.0 /
           static_cast<double>(freq.QuadPart);
#else
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec * 1000.0 + ts.tv_nsec / 1000000.0;
#endif
}

static long get_file_size(const char* filename) {
    struct stat st;
    if (stat(filename, &st) == 0) {
        return static_cast<long>(st.st_size);
    }
    return 0;
}

static const char* get_temp_dir(void) {
    const char* override = std::getenv("CARQUET_BENCH_TMPDIR");
    if (override && override[0] != '\0') {
        return override;
    }
#ifdef _WIN32
    static char temp_dir[512] = {0};
    if (temp_dir[0] == 0) {
        const char* tmp = std::getenv("TEMP");
        if (!tmp) tmp = std::getenv("TMP");
        if (!tmp) tmp = ".";
        std::snprintf(temp_dir, sizeof(temp_dir), "%s", tmp);
    }
    return temp_dir;
#else
    return "/tmp";
#endif
}

static void purge_file_cache(const char* filename) {
#if defined(__APPLE__)
    char tmp[512];
    std::snprintf(tmp, sizeof(tmp), "%s.nocache", filename);

    int src = open(filename, O_RDONLY);
    if (src < 0) return;

    int dst = open(tmp, O_WRONLY | O_CREAT | O_TRUNC, 0644);
    if (dst < 0) {
        close(src);
        return;
    }

    fcntl(dst, F_NOCACHE, 1);
    fcntl(src, F_NOCACHE, 1);

    char buf[262144];
    ssize_t n = 0;
    while ((n = read(src, buf, sizeof(buf))) > 0) {
        (void)write(dst, buf, static_cast<size_t>(n));
    }

    close(src);
    close(dst);
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

static uint32_t lcg_state = 42;

static void lcg_seed(uint32_t seed) {
    lcg_state = seed;
}

static uint32_t lcg_rand(void) {
    lcg_state = lcg_state * 1103515245u + 12345u;
    return (lcg_state >> 16) & 0x7FFFu;
}

static double lcg_normal(void) {
    double u1 = (lcg_rand() + 1.0) / 32768.0;
    double u2 = (lcg_rand() + 1.0) / 32768.0;
    return std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * 3.14159265358979 * u2);
}

static bool check_status(const char* what, const arrow::Status& status) {
    if (status.ok()) {
        return true;
    }
    std::fprintf(stderr, "%s: %s\n", what, status.ToString().c_str());
    return false;
}

template <typename T>
static bool take_result(const char* what, arrow::Result<T> result, T* out) {
    if (!result.ok()) {
        std::fprintf(stderr, "%s: %s\n", what, result.status().ToString().c_str());
        return false;
    }
    *out = result.MoveValueUnsafe();
    return true;
}

static bool build_file_reader(const char* filename, bool use_threads,
                              std::unique_ptr<parquet::arrow::FileReader>* out) {
    parquet::ReaderProperties reader_props;
    reader_props.set_page_checksum_verification(false);

    parquet::ArrowReaderProperties arrow_props;
    arrow_props.set_use_threads(use_threads);
    arrow_props.set_pre_buffer(true);

    parquet::arrow::FileReaderBuilder builder;
    if (!check_status("Open input file",
                      builder.OpenFile(filename, true, reader_props))) {
        return false;
    }
    builder.properties(arrow_props);

    if (!take_result("Build file reader", builder.Build(), out)) {
        return false;
    }

    (*out)->set_use_threads(use_threads);
    return true;
}

static std::unique_ptr<test_data_t> test_data_create(int num_rows, bool nullable) {
    std::unique_ptr<arrow::Buffer> ids_buf_tmp;
    std::unique_ptr<arrow::Buffer> values_buf_tmp;
    std::unique_ptr<arrow::Buffer> categories_buf_tmp;

    if (!take_result("Allocate ids buffer",
                     arrow::AllocateBuffer(num_rows * static_cast<int64_t>(sizeof(int64_t))),
                     &ids_buf_tmp)) {
        return nullptr;
    }
    if (!take_result("Allocate values buffer",
                     arrow::AllocateBuffer(num_rows * static_cast<int64_t>(sizeof(double))),
                     &values_buf_tmp)) {
        return nullptr;
    }
    if (!take_result("Allocate categories buffer",
                     arrow::AllocateBuffer(num_rows * static_cast<int64_t>(sizeof(int32_t))),
                     &categories_buf_tmp)) {
        return nullptr;
    }

    std::shared_ptr<arrow::Buffer> ids_buf(std::move(ids_buf_tmp));
    std::shared_ptr<arrow::Buffer> values_buf(std::move(values_buf_tmp));
    std::shared_ptr<arrow::Buffer> categories_buf(std::move(categories_buf_tmp));

    int64_t* ids = reinterpret_cast<int64_t*>(ids_buf->mutable_data());
    double* values = reinterpret_cast<double*>(values_buf->mutable_data());
    int32_t* categories = reinterpret_cast<int32_t*>(categories_buf->mutable_data());

    lcg_seed(42);
    for (int i = 0; i < num_rows; ++i) {
        ids[i] = 1000000 + static_cast<int64_t>(lcg_rand() % 9000000u);
        values[i] = std::fabs(100.0 + 50.0 * lcg_normal());
        categories[i] = static_cast<int32_t>(lcg_rand() % 100u);
    }

    auto td = std::make_unique<test_data_t>();
    td->num_rows = num_rows;
    td->nullable = nullable;

    /* Nullable datasets: same values plus an Arrow validity bitmap per column
     * (Arrow's native input; the values under null slots are ignored). */
    std::shared_ptr<arrow::Buffer> validity[NUM_COLUMNS];
    for (int c = 0; c < NUM_COLUMNS; ++c) {
        td->null_counts[c] = 0;
        if (!nullable) {
            continue;
        }
        std::shared_ptr<arrow::Buffer> bitmap;
        if (!take_result("Allocate validity bitmap", arrow::AllocateBitmap(num_rows),
                         &bitmap)) {
            return nullptr;
        }
        uint8_t* bits = bitmap->mutable_data();
        std::memset(bits, 0, static_cast<size_t>(bitmap->size()));
        for (int i = 0; i < num_rows; ++i) {
            if (bench_is_null(static_cast<uint32_t>(i), c)) {
                td->null_counts[c]++;
            } else {
                bits[i >> 3] = static_cast<uint8_t>(bits[i >> 3] | (1u << (i & 7)));
            }
        }
        validity[c] = std::move(bitmap);
    }

    auto ids_array = std::make_shared<arrow::Int64Array>(
        num_rows, ids_buf, validity[0], td->null_counts[0]);
    auto values_array = std::make_shared<arrow::DoubleArray>(
        num_rows, values_buf, validity[1], td->null_counts[1]);
    auto categories_array = std::make_shared<arrow::Int32Array>(
        num_rows, categories_buf, validity[2], td->null_counts[2]);

    td->table = arrow::Table::Make(
        arrow::schema({
            arrow::field("id", arrow::int64(), nullable),
            arrow::field("value", arrow::float64(), nullable),
            arrow::field("category", arrow::int32(), nullable),
        }),
        {ids_array, values_array, categories_array},
        num_rows);
    return td;
}

static double trimmed_median(const double* values, int n) {
    std::vector<double> sorted(values, values + n);
    std::sort(sorted.begin(), sorted.end());
    if (n <= 2) {
        return sorted[n / 2];
    }
    int lo = 1;
    int hi = n - 2;
    int trimmed = hi - lo + 1;
    if ((trimmed & 1) != 0) {
        return sorted[lo + trimmed / 2];
    }
    return (sorted[lo + trimmed / 2 - 1] + sorted[lo + trimmed / 2]) / 2.0;
}

template <typename ArrayType>
static void accumulate_first_last(const std::shared_ptr<arrow::ChunkedArray>& column,
                                  volatile int64_t* checksum) {
    for (const auto& chunk : column->chunks()) {
        auto array = std::static_pointer_cast<ArrayType>(chunk);
        if (array->length() <= 0) {
            continue;
        }
        *checksum += static_cast<int64_t>(array->Value(0));
        *checksum += static_cast<int64_t>(array->Value(array->length() - 1));
    }
}

static double benchmark_write(const char* filename, const test_data_t* td,
                              parquet::Compression::type codec) {
    double start = get_time_ms();

    std::shared_ptr<arrow::io::FileOutputStream> sink;
    if (!take_result("Open output file",
                     arrow::io::FileOutputStream::Open(filename, false), &sink)) {
        return -1.0;
    }

    parquet::WriterProperties::Builder props_builder;
    props_builder.disable_dictionary();
    if (td->nullable) {
        /* Nullable datasets: dictionary on `category` only, like carquet's
         * carquet_writer_set_column_encoding(w, 2, RLE_DICTIONARY). */
        props_builder.enable_dictionary("category");
    }
    props_builder.enable_page_checksum();
    props_builder.compression(codec);

    if (codec != parquet::Compression::UNCOMPRESSED) {
        props_builder.encoding("value", parquet::Encoding::BYTE_STREAM_SPLIT);
    }
    if (codec == parquet::Compression::ZSTD) {
        props_builder.compression_level(get_benchmark_zstd_level());
    }

    parquet::ArrowWriterProperties::Builder arrow_props_builder;
    auto props = props_builder.build();
    auto arrow_props = arrow_props_builder.build();

    int64_t rg_rows = td->num_rows > 1000000 ? td->num_rows / 10 : 100000;

    arrow::Status status = parquet::arrow::WriteTable(
        *td->table,
        arrow::default_memory_pool(),
        sink,
        rg_rows,
        props,
        arrow_props);
    if (!check_status("WriteTable", status)) {
        (void)sink->Close();
        return -1.0;
    }
    if (!check_status("Close output file", sink->Close())) {
        return -1.0;
    }
    return get_time_ms() - start;
}

/**
 * Low-level Parquet read benchmark.
 *
 * Uses parquet::ParquetFileReader directly (NOT the Arrow conversion layer).
 * This measures pure Parquet I/O + decompression + decoding — the same work
 * that Carquet's batch reader does — without the overhead of materializing
 * Arrow Tables, allocating Arrow arrays, or converting between type systems.
 *
 * Row groups are read in parallel using std::async, matching the parallelism
 * strategy that Carquet's pipeline uses.
 */

struct rg_read_result_t {
    int64_t rows;
    int64_t checksum;
    int64_t nulls[NUM_COLUMNS];
};

/* OPTIONAL columns (max_definition_level > 0) are read with a def-level
 * buffer: ReadBatch then returns dense non-null values plus one level per
 * row, and the nulls are counted from the levels (Arrow's cheapest
 * low-level form; carquet additionally spaces values and builds a bitmap).
 * REQUIRED columns take the exact same path as before. */
static rg_read_result_t read_row_group_lowlevel(
    const char* filename, int rg_index, int64_t rg_rows, int num_columns,
    bool use_mmap) {

    rg_read_result_t result = {0, 0, {0, 0, 0}};

    parquet::ReaderProperties props = parquet::default_reader_properties();
    props.set_page_checksum_verification(false);
    props.enable_buffered_stream();
    props.set_buffer_size(4 * 1024 * 1024);

    std::unique_ptr<parquet::ParquetFileReader> reader;
    try {
        reader = parquet::ParquetFileReader::OpenFile(filename, use_mmap, props);
    } catch (...) {
        return result;
    }

    auto rg_reader = reader->RowGroup(rg_index);

    std::vector<int64_t> col0_buf(rg_rows);
    std::vector<double>  col1_buf(rg_rows);
    std::vector<int32_t> col2_buf(rg_rows);
    std::vector<int16_t> def_buf;
    int64_t values_total[NUM_COLUMNS] = {0, 0, 0};

    for (int col = 0; col < num_columns && col < 3; col++) {
        auto col_reader = rg_reader->Column(col);
        int16_t max_def = col_reader->descr()->max_definition_level();
        int16_t* def_levels = nullptr;
        if (max_def > 0) {
            def_buf.resize(static_cast<size_t>(rg_rows));
            def_levels = def_buf.data();
        }
        int64_t values_read = 0;
        int64_t rows_remaining = rg_rows;
        int64_t values_offset = 0;

        while (col_reader->HasNext() && rows_remaining > 0) {
            int64_t batch;
            int16_t* defs = def_levels ? def_levels + (rg_rows - rows_remaining) : nullptr;
            if (col == 0) {
                batch = static_cast<parquet::Int64Reader*>(col_reader.get())
                    ->ReadBatch(rows_remaining, defs, nullptr,
                                col0_buf.data() + values_offset,
                                &values_read);
            } else if (col == 1) {
                batch = static_cast<parquet::DoubleReader*>(col_reader.get())
                    ->ReadBatch(rows_remaining, defs, nullptr,
                                col1_buf.data() + values_offset,
                                &values_read);
            } else {
                batch = static_cast<parquet::Int32Reader*>(col_reader.get())
                    ->ReadBatch(rows_remaining, defs, nullptr,
                                col2_buf.data() + values_offset,
                                &values_read);
            }
            rows_remaining -= batch;
            values_offset += values_read;
        }
        values_total[col] = values_offset;

        if (def_levels) {
            int64_t rows_read = rg_rows - rows_remaining;
            int64_t nulls = 0;
            for (int64_t i = 0; i < rows_read; i++) {
                nulls += def_levels[i] < max_def;
            }
            result.nulls[col] = nulls;
        }
    }

    if (values_total[0] > 0) result.checksum += col0_buf[0] + col0_buf[values_total[0] - 1];
    if (values_total[1] > 0) {
        result.checksum += static_cast<int64_t>(col1_buf[0] + col1_buf[values_total[1] - 1]);
    }
    if (values_total[2] > 0) result.checksum += col2_buf[0] + col2_buf[values_total[2] - 1];
    result.rows = rg_rows;
    return result;
}

/* Returns the elapsed ms, or -1.0 on failure (including a null count that
 * does not match expected_nulls when it is given). */
static double benchmark_read(const char* filename, int expected_rows, bool use_mmap,
                             const int64_t* expected_nulls) {
    double start = get_time_ms();

    /* Open once to read metadata */
    parquet::ReaderProperties props = parquet::default_reader_properties();
    props.set_page_checksum_verification(false);

    std::unique_ptr<parquet::ParquetFileReader> reader;
    try {
        reader = parquet::ParquetFileReader::OpenFile(filename, use_mmap, props);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "Failed to open %s: %s\n", filename, e.what());
        return -1.0;
    }

    auto file_meta = reader->metadata();
    int num_row_groups = file_meta->num_row_groups();
    int num_columns = file_meta->num_columns();

    /* Collect row group sizes before closing the metadata reader */
    std::vector<int64_t> rg_sizes(num_row_groups);
    for (int i = 0; i < num_row_groups; i++)
        rg_sizes[i] = file_meta->RowGroup(i)->num_rows();
    reader.reset();

    /* Read row groups in parallel — each opens its own file reader
     * (independent mmap/buffered I/O, matching Carquet's per-slot mmap). */
    int64_t total_rows = 0;
    int64_t total_nulls[NUM_COLUMNS] = {0, 0, 0};
    volatile int64_t checksum = 0;

    if (num_row_groups > 1) {
        unsigned hw = std::thread::hardware_concurrency();
        if (hw == 0) hw = 4;
        int max_parallel = static_cast<int>(hw);
        if (max_parallel > num_row_groups) max_parallel = num_row_groups;

        std::vector<std::future<rg_read_result_t>> futures;
        futures.reserve(num_row_groups);
        for (int rg = 0; rg < num_row_groups; rg++) {
            futures.push_back(std::async(
                std::launch::async,
                read_row_group_lowlevel, filename, rg, rg_sizes[rg], num_columns,
                use_mmap));
        }
        for (auto& f : futures) {
            auto r = f.get();
            total_rows += r.rows;
            checksum += r.checksum;
            for (int c = 0; c < NUM_COLUMNS; c++) total_nulls[c] += r.nulls[c];
        }
    } else {
        /* Single row group — read directly, no async overhead */
        auto r = read_row_group_lowlevel(filename, 0, rg_sizes[0], num_columns, use_mmap);
        total_rows = r.rows;
        checksum = r.checksum;
        for (int c = 0; c < NUM_COLUMNS; c++) total_nulls[c] = r.nulls[c];
    }
    (void)checksum;

    if (total_rows != expected_rows) {
        std::fprintf(stderr, "Warning: row count mismatch %lld vs %d\n",
                     static_cast<long long>(total_rows), expected_rows);
    }

    double elapsed = get_time_ms() - start;
    if (expected_nulls) {
        for (int c = 0; c < NUM_COLUMNS; c++) {
            if (total_nulls[c] != expected_nulls[c]) {
                std::fprintf(stderr, "Null count mismatch in column %d of %s: %lld vs %lld\n",
                             c, filename, static_cast<long long>(total_nulls[c]),
                             static_cast<long long>(expected_nulls[c]));
                return -1.0;
            }
        }
    }
    return elapsed;
}

static bool run_benchmark(const char* dataset_name, int num_rows, bool nullable,
                          parquet::Compression::type codec,
                          const char* compression_name) {
    char filename[256];
    std::snprintf(filename, sizeof(filename), "%s/benchmark_%s_%s_arrow_cpp.parquet",
                  get_temp_dir(), dataset_name, compression_name);

    std::printf("\n=== %s (%d rows, %s) ===\n", dataset_name, num_rows, compression_name);

    std::unique_ptr<test_data_t> td = test_data_create(num_rows, nullable);
    if (!td) {
        std::fprintf(stderr, "Failed to create test data\n");
        return false;
    }

    int iters;
    const char* iter_env = std::getenv("CARQUET_BENCH_ITERATIONS");
    if (iter_env && std::atoi(iter_env) > 0) {
        iters = std::atoi(iter_env);
        if (iters > MAX_BENCH_ITERATIONS) iters = MAX_BENCH_ITERATIONS;
    } else if (num_rows <= 100000) {
        iters = BENCH_ITERATIONS_SMALL;
    } else if (num_rows <= 1000000) {
        iters = BENCH_ITERATIONS_MEDIUM;
    } else {
        iters = BENCH_ITERATIONS_LARGE;
    }

    for (int i = 0; i < WARMUP_ITERATIONS; ++i) {
        if (benchmark_write(filename, td.get(), codec) <= 0.0) {
            std::remove(filename);
            return false;
        }
    }

    double write_times[MAX_BENCH_ITERATIONS];
    for (int i = 0; i < iters; ++i) {
        write_times[i] = benchmark_write(filename, td.get(), codec);
        if (write_times[i] <= 0.0) {
            std::remove(filename);
            return false;
        }
    }

    long file_size = get_file_size(filename);
    if (file_size <= 0) {
        std::fprintf(stderr, "Benchmark file is empty or missing\n");
        std::remove(filename);
        return false;
    }
    const int64_t* expected_nulls = nullable ? td->null_counts : nullptr;
    purge_file_cache(filename);

    for (int i = 0; i < WARMUP_ITERATIONS; ++i) {
        if (benchmark_read(filename, num_rows, true, expected_nulls) <= 0.0) {
            std::remove(filename);
            return false;
        }
    }

    double read_times[MAX_BENCH_ITERATIONS];
    for (int i = 0; i < iters; ++i) {
        read_times[i] = benchmark_read(filename, num_rows, true, expected_nulls);
        if (read_times[i] <= 0.0) {
            std::remove(filename);
            return false;
        }
    }

    /* Same protocol without memory mapping (Arrow's default), after the mmap
     * pass so that one is unchanged. */
    purge_file_cache(filename);

    for (int i = 0; i < WARMUP_ITERATIONS; ++i) {
        if (benchmark_read(filename, num_rows, false, expected_nulls) <= 0.0) {
            std::remove(filename);
            return false;
        }
    }

    double read_nommap_times[MAX_BENCH_ITERATIONS];
    for (int i = 0; i < iters; ++i) {
        read_nommap_times[i] = benchmark_read(filename, num_rows, false, expected_nulls);
        if (read_nommap_times[i] <= 0.0) {
            std::remove(filename);
            return false;
        }
    }

    double write_med = trimmed_median(write_times, iters);
    double read_med = trimmed_median(read_times, iters);
    double read_nommap_med = trimmed_median(read_nommap_times, iters);
    double rows_per_sec_write = (num_rows / write_med) * 1000.0;
    double rows_per_sec_read = (num_rows / read_med) * 1000.0;

    std::printf("  Write: %.2f ms (%.2f M rows/sec)\n",
                write_med, rows_per_sec_write / 1e6);
    std::printf("  Read:  %.2f ms (%.2f M rows/sec)\n",
                read_med, rows_per_sec_read / 1e6);
    std::printf("  Read (no mmap): %.2f ms (%.2f M rows/sec)\n",
                read_nommap_med, (num_rows / read_nommap_med) * 1000.0 / 1e6);
    std::printf("  File:  %.2f MB (%.2f bytes/row)\n",
                file_size / (1024.0 * 1024.0), static_cast<double>(file_size) / num_rows);

    /* 8th field (read_nommap) is optional for the orchestrator's parser. */
    std::printf("CSV:arrow_cpp,%s,%s,%d,%.2f,%.2f,%ld,%.2f\n",
                dataset_name, compression_name, num_rows, write_med, read_med, file_size,
                read_nommap_med);

    /* File is NOT removed here — the orchestrator handles cleanup
     * after cross-read benchmarking is done. */
    return true;
}

static int find_dataset(const char* name, dataset_t* datasets, int n) {
    for (int i = 0; i < n; ++i) {
        if (std::strcmp(datasets[i].name, name) == 0) {
            return i;
        }
    }
    return -1;
}

static int find_compression(const char* name, compression_config_t* comps, int n) {
    for (int i = 0; i < n; ++i) {
        if (std::strcmp(comps[i].name, name) == 0) {
            return i;
        }
    }
    return -1;
}

static int run_crossread(const char* filename) {
    /* Open file to get row count from metadata */
    std::unique_ptr<parquet::arrow::FileReader> reader;
    if (!build_file_reader(filename, true, &reader)) {
        std::fprintf(stderr, "Failed to open %s\n", filename);
        return 1;
    }

    int64_t num_rows = reader->parquet_reader()->metadata()->num_rows();
    reader.reset();

    if (num_rows <= 0) {
        std::fprintf(stderr, "File %s has no rows\n", filename);
        return 1;
    }

    int expected_rows = static_cast<int>(num_rows);
    long file_size = get_file_size(filename);

    /* Scale iterations based on row count */
    int iters;
    const char* iter_env = std::getenv("CARQUET_BENCH_ITERATIONS");
    if (iter_env && std::atoi(iter_env) > 0) {
        iters = std::atoi(iter_env);
        if (iters > MAX_BENCH_ITERATIONS) iters = MAX_BENCH_ITERATIONS;
    } else if (num_rows <= 100000) {
        iters = BENCH_ITERATIONS_SMALL;
    } else if (num_rows <= 1000000) {
        iters = BENCH_ITERATIONS_MEDIUM;
    } else {
        iters = BENCH_ITERATIONS_LARGE;
    }

    /* Purge file cache before reads */
    purge_file_cache(filename);

    /* Warmup reads */
    for (int i = 0; i < WARMUP_ITERATIONS; ++i) {
        double t = benchmark_read(filename, expected_rows, true, nullptr);
        if (t <= 0.0) {
            std::fprintf(stderr, "Warmup read failed for %s\n", filename);
            return 1;
        }
    }

    /* Measured reads */
    double read_times[MAX_BENCH_ITERATIONS];
    for (int i = 0; i < iters; ++i) {
        read_times[i] = benchmark_read(filename, expected_rows, true, nullptr);
        if (read_times[i] <= 0.0) {
            std::fprintf(stderr, "Benchmark read failed for %s\n", filename);
            return 1;
        }
    }

    double read_med = trimmed_median(read_times, iters);

    std::printf("CSV:arrow_cpp,crossread,external,%d,0.00,%.2f,%ld\n",
                expected_rows, read_med, file_size);

    return 0;
}

int main(int argc, char* argv[]) {
    std::setvbuf(stdout, NULL, _IONBF, 0);

    if (argc == 2 && std::strcmp(argv[1], "--version") == 0) {
        std::printf("%s\n", ARROW_VERSION_STRING);
        return 0;
    }

    /* read <filename> subcommand for cross-read benchmarking */
    if (argc == 2 && std::strcmp(argv[1], "read") == 0) {
        std::fprintf(stderr, "Usage: %s read <filename>\n", argv[0]);
        return 1;
    }
    if (argc >= 3 && std::strcmp(argv[1], "read") == 0) {
        return run_crossread(argv[2]);
    }

    dataset_t datasets[] = {
        {"small", 100000, false},
        {"medium", 1000000, false},
        {"large", 10000000, false},
        {"xlarge", 100000000, false},
        /* OPTIONAL columns, ~NULL_PERCENT % nulls, dictionary on category */
        {"small_nullable", 100000, true},
        {"medium_nullable", 1000000, true},
        {"large_nullable", 10000000, true},
        {"xlarge_nullable", 100000000, true}
    };
    int num_datasets = sizeof(datasets) / sizeof(datasets[0]);

    compression_config_t compressions[] = {
        {parquet::Compression::UNCOMPRESSED, "none"},
        {parquet::Compression::SNAPPY, "snappy"},
        {parquet::Compression::ZSTD, "zstd"},
        {parquet::Compression::LZ4, "lz4"}
    };
    int num_compressions = sizeof(compressions) / sizeof(compressions[0]);

    if (argc == 3) {
        int d = find_dataset(argv[1], datasets, num_datasets);
        int c = find_compression(argv[2], compressions, num_compressions);
        if (d < 0 || c < 0) {
            std::fprintf(stderr, "Usage: %s [dataset] [compression]\n", argv[0]);
            std::fprintf(stderr, "  datasets:     small, medium, large, xlarge (+ _nullable)\n");
            std::fprintf(stderr, "  compressions: none, snappy, zstd, lz4\n");
            return 1;
        }
        return run_benchmark(datasets[d].name, datasets[d].rows, datasets[d].nullable,
                             compressions[c].codec, compressions[c].name) ? 0 : 1;
    }

    std::printf("Arrow C++ Benchmark\n");

    for (int d = 0; d < num_datasets; ++d) {
        if (datasets[d].nullable) continue;  /* opt-in: pass "<size>_nullable <codec>" */
        for (int c = 0; c < num_compressions; ++c) {
            if (!run_benchmark(datasets[d].name, datasets[d].rows, datasets[d].nullable,
                               compressions[c].codec, compressions[c].name)) {
                return 1;
            }
        }
    }

    std::printf("\nBenchmark complete.\n");
    return 0;
}
