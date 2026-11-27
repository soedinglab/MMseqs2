/**
 * @file compression_codecs.c
 * @brief Example demonstrating different compression codecs
 *
 * This example demonstrates:
 * - Writing the same data with different compression codecs
 * - Comparing file sizes across codecs
 * - Reading compressed data back
 *
 * Supported codecs:
 * - UNCOMPRESSED
 * - SNAPPY
 * - GZIP
 * - LZ4_RAW (the LZ4 codec to use for new files)
 * - LZ4 (deprecated Hadoop-framed variant, kept for compatibility)
 * - ZSTD
 * - BROTLI (when built with CARQUET_ENABLE_BROTLI, the default)
 *
 * Build:
 *   gcc -o compression_codecs compression_codecs.c -lcarquet -I../include
 *
 * Run:
 *   ./compression_codecs
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include <carquet/carquet.h>
#include "example_paths.h"

#define NUM_ROWS 10000

/**
 * Get file size
 */
static int64_t get_file_size(const char* filename) {
    struct stat st;
    if (stat(filename, &st) == 0) {
        return st.st_size;
    }
    return -1;
}

/**
 * Create a simple schema
 */
static carquet_schema_t* create_schema(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema) {
        fprintf(stderr, "Failed to create schema\n");
        return NULL;
    }

    carquet_logical_type_t string_type = { .id = CARQUET_LOGICAL_STRING };

(void)carquet_schema_add_column(schema, "id",
        CARQUET_PHYSICAL_INT64, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
(void)carquet_schema_add_column(schema, "value",
        CARQUET_PHYSICAL_DOUBLE, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
(void)carquet_schema_add_column(schema, "category",
        CARQUET_PHYSICAL_BYTE_ARRAY, &string_type, CARQUET_REPETITION_REQUIRED, 0, 0);

    return schema;
}

/**
 * Generate test data
 */
typedef struct {
    int64_t* ids;
    double* values;
    carquet_byte_array_t* categories;
    char category_buffers[NUM_ROWS][32];
} test_data_t;

static test_data_t* generate_test_data(void) {
    test_data_t* data = malloc(sizeof(test_data_t));
    if (!data) return NULL;

    data->ids = malloc(NUM_ROWS * sizeof(int64_t));
    data->values = malloc(NUM_ROWS * sizeof(double));
    data->categories = malloc(NUM_ROWS * sizeof(carquet_byte_array_t));

    if (!data->ids || !data->values || !data->categories) {
        free(data->ids);
        free(data->values);
        free(data->categories);
        free(data);
        return NULL;
    }

    /* Categories that repeat (good for dictionary encoding) */
    const char* category_names[] = {
        "electronics", "clothing", "food", "furniture", "books",
        "sports", "toys", "health", "automotive", "garden"
    };
    const int num_categories = 10;

    for (int i = 0; i < NUM_ROWS; i++) {
        data->ids[i] = (int64_t)i + 1;

        /* Values with patterns (good for compression) */
        data->values[i] = 100.0 + (i % 100) * 1.5 + (i / 1000) * 10.0;

        /* Repeating categories */
        const char* cat = category_names[i % num_categories];
        strcpy(data->category_buffers[i], cat);
        data->categories[i].data = (uint8_t*)data->category_buffers[i];
        data->categories[i].length = (int32_t)strlen(cat);
    }

    return data;
}

static void free_test_data(test_data_t* data) {
    if (data) {
        free(data->ids);
        free(data->values);
        free(data->categories);
        free(data);
    }
}

/**
 * Write file with specified compression
 */
static int write_with_compression(
    const char* filename,
    carquet_schema_t* schema,
    test_data_t* data,
    carquet_compression_t codec)
{
    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = codec;

    carquet_writer_t* writer = carquet_writer_create(filename, schema, &opts, &err);
    if (!writer) {
        fprintf(stderr, "    Failed to create writer: %s\n", err.message);
        return -1;
    }

    carquet_status_t status;

    status = carquet_writer_write_batch(writer, 0, data->ids, NUM_ROWS, NULL, NULL);
    if (status != CARQUET_OK) goto error;

    status = carquet_writer_write_batch(writer, 1, data->values, NUM_ROWS, NULL, NULL);
    if (status != CARQUET_OK) goto error;

    status = carquet_writer_write_batch(writer, 2, data->categories, NUM_ROWS, NULL, NULL);
    if (status != CARQUET_OK) goto error;

    status = carquet_writer_close(writer);
    if (status != CARQUET_OK) {
        fprintf(stderr, "    Failed to close writer\n");
        return -1;
    }

    return 0;

error:
    fprintf(stderr, "    Failed to write data\n");
    carquet_writer_abort(writer);
    return -1;
}

/**
 * Verify file can be read
 */
static int verify_file(const char* filename) {
    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_reader_t* reader = carquet_reader_open(filename, NULL, &err);
    if (!reader) {
        fprintf(stderr, "    Failed to open for verification: %s\n", err.message);
        return -1;
    }

    int64_t num_rows = carquet_reader_num_rows(reader);
    if (num_rows != NUM_ROWS) {
        fprintf(stderr, "    Row count mismatch: expected %d, got %lld\n",
                NUM_ROWS, (long long)num_rows);
        carquet_reader_close(reader);
        return -1;
    }

    /* Read first few values to verify */
    carquet_column_reader_t* col = carquet_reader_get_column(reader, 0, 0, &err);
    if (col) {
        int64_t ids[10];
        int64_t count = carquet_column_read_batch(col, ids, 10, NULL, NULL);
        if (count > 0 && ids[0] != 1) {
            fprintf(stderr, "    Data verification failed\n");
            carquet_column_reader_free(col);
            carquet_reader_close(reader);
            return -1;
        }
        carquet_column_reader_free(col);
    }

    carquet_reader_close(reader);
    return 0;
}

/**
 * Test a single compression codec
 */
static int test_codec(
    carquet_compression_t codec,
    const char* base_path,
    carquet_schema_t* schema,
    test_data_t* data,
    int64_t* out_size)
{
    const char* codec_name = carquet_compression_name(codec);

    char filename[256];
    snprintf(filename, sizeof(filename), "%s_%s.parquet",
             base_path, codec_name);

    printf("  Testing %s...\n", codec_name);

    /* Write file */
    if (write_with_compression(filename, schema, data, codec) != 0) {
        printf("    FAILED to write\n");
        *out_size = -1;
        remove(filename);
        return 1;
    }

    /* Get file size */
    int64_t size = get_file_size(filename);
    *out_size = size;

    /* Verify file */
    if (verify_file(filename) != 0) {
        printf("    FAILED verification\n");
        remove(filename);
        return 1;
    }

    printf("    Size: %lld bytes (%.2f KB)\n",
           (long long)size, size / 1024.0);

    /* Clean up */
    remove(filename);
    return 0;
}

int main(int argc, char* argv[]) {
    printf("=== Carquet Compression Codecs Example ===\n");
    printf("Library version: %s\n\n", carquet_version());

    /* Initialize */
(void)carquet_init();

    /* Create schema */
    carquet_schema_t* schema = create_schema();
    if (!schema) {
        return 1;
    }

    /* Generate test data */
    printf("Generating %d rows of test data...\n", NUM_ROWS);
    test_data_t* data = generate_test_data();
    if (!data) {
        carquet_schema_free(schema);
        fprintf(stderr, "Failed to generate test data\n");
        return 1;
    }

    /* Calculate raw data size */
    int64_t raw_size = NUM_ROWS * sizeof(int64_t)      /* ids */
                     + NUM_ROWS * sizeof(double)        /* values */
                     + NUM_ROWS * 12;                   /* average category length */
    printf("Raw data size: ~%lld bytes (%.2f KB)\n\n", (long long)raw_size, raw_size / 1024.0);

    /* Test each codec */
    printf("Writing files with different compression codecs:\n\n");

    char base_path[1024];
    example_temp_path(base_path, sizeof(base_path), "compression_test");

    typedef struct {
        carquet_compression_t codec;
        int64_t size;
    } codec_result_t;

    codec_result_t results[] = {
        { CARQUET_COMPRESSION_UNCOMPRESSED, 0 },
        { CARQUET_COMPRESSION_SNAPPY, 0 },
        { CARQUET_COMPRESSION_GZIP, 0 },
        { CARQUET_COMPRESSION_LZ4_RAW, 0 },
        { CARQUET_COMPRESSION_LZ4, 0 },
        { CARQUET_COMPRESSION_ZSTD, 0 },
        { CARQUET_COMPRESSION_BROTLI, 0 },
    };
    const int num_codecs = sizeof(results) / sizeof(results[0]);

    int failures = 0;
    for (int i = 0; i < num_codecs; i++) {
        int rc = test_codec(results[i].codec, base_path, schema, data, &results[i].size);
        /* BROTLI is optional at build time (CARQUET_ENABLE_BROTLI=OFF). */
        if (rc != 0 && results[i].codec == CARQUET_COMPRESSION_BROTLI) {
            printf("    (BROTLI not available in this build; skipped)\n");
        } else {
            failures += rc;
        }
    }

    /* Summary */
    printf("\n=== Compression Summary ===\n\n");
    printf("%-15s %12s %12s\n", "Codec", "Size (KB)", "Ratio");
    printf("%-15s %12s %12s\n", "-----", "---------", "-----");

    int64_t uncompressed_size = results[0].size;
    for (int i = 0; i < num_codecs; i++) {
        if (results[i].size > 0) {
            double ratio = 100.0 * results[i].size / uncompressed_size;
            printf("%-15s %12.2f %11.1f%%\n",
                   carquet_compression_name(results[i].codec),
                   results[i].size / 1024.0,
                   ratio);
        }
    }

    /* Recommendations */
    printf("\nCodec Recommendations:\n");
    printf("  - SNAPPY: Fast compression/decompression, good for real-time workloads\n");
    printf("  - ZSTD: Best compression ratio, good balance of speed and size\n");
    printf("  - LZ4_RAW: Fastest decompression, good for read-heavy workloads\n");
    printf("    (LZ4, the Hadoop-framed codec id, is deprecated: write LZ4_RAW)\n");
    printf("  - BROTLI: High ratio on text-like data, slower to compress\n");
    printf("  - GZIP: Widely compatible, good compression but slower\n");

    /* Clean up */
    free_test_data(data);
    carquet_schema_free(schema);

    if (failures) {
        printf("\n=== Compression example: %d codec(s) FAILED ===\n", failures);
        return 1;
    }
    printf("\n=== Compression example completed ===\n");
    return 0;
}
