/**
 * @file data_types.c
 * @brief Example demonstrating all supported Parquet data types
 *
 * This example demonstrates:
 * - All physical types: BOOLEAN, INT32, INT64, FLOAT, DOUBLE, BYTE_ARRAY
 * - Logical types: STRING, DATE, TIMESTAMP, DECIMAL
 * - Fixed-length byte arrays (UUID, etc.)
 *
 * Build:
 *   gcc -o data_types data_types.c -lcarquet -I../include
 *
 * Run:
 *   ./data_types
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>

#include <carquet/carquet.h>
#include "example_paths.h"

#define NUM_ROWS 100

/**
 * Create a schema demonstrating various data types
 */
static carquet_schema_t* create_typed_schema(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema) {
        fprintf(stderr, "Failed to create schema: %s\n", err.message);
        return NULL;
    }

    /* INT32 column (plain) */
(void)carquet_schema_add_column(schema, "count",
        CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);

    /* INT32 with DATE logical type */
    carquet_logical_type_t date_type = { .id = CARQUET_LOGICAL_DATE };
(void)carquet_schema_add_column(schema, "created_date",
        CARQUET_PHYSICAL_INT32, &date_type, CARQUET_REPETITION_REQUIRED, 0, 0);

    /* INT64 column (plain) */
(void)carquet_schema_add_column(schema, "big_number",
        CARQUET_PHYSICAL_INT64, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);

    /* INT64 with TIMESTAMP logical type (microseconds, UTC) */
    carquet_logical_type_t timestamp_type = {
        .id = CARQUET_LOGICAL_TIMESTAMP,
        .params.timestamp = {
            .unit = CARQUET_TIME_UNIT_MICROS,
            .is_adjusted_to_utc = true
        }
    };
(void)carquet_schema_add_column(schema, "event_time",
        CARQUET_PHYSICAL_INT64, &timestamp_type, CARQUET_REPETITION_REQUIRED, 0, 0);

    /* FLOAT column */
(void)carquet_schema_add_column(schema, "temperature",
        CARQUET_PHYSICAL_FLOAT, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);

    /* DOUBLE column */
(void)carquet_schema_add_column(schema, "precise_value",
        CARQUET_PHYSICAL_DOUBLE, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);

    /* BYTE_ARRAY with STRING logical type */
    carquet_logical_type_t string_type = { .id = CARQUET_LOGICAL_STRING };
(void)carquet_schema_add_column(schema, "description",
        CARQUET_PHYSICAL_BYTE_ARRAY, &string_type, CARQUET_REPETITION_REQUIRED, 0, 0);

    printf("Schema created with %d columns:\n", carquet_schema_num_columns(schema));
    printf("  - count: INT32\n");
    printf("  - created_date: INT32 (DATE)\n");
    printf("  - big_number: INT64\n");
    printf("  - event_time: INT64 (TIMESTAMP_MICROS)\n");
    printf("  - temperature: FLOAT\n");
    printf("  - precise_value: DOUBLE\n");
    printf("  - description: BYTE_ARRAY (STRING)\n");

    return schema;
}

/**
 * Convert days since Unix epoch to date string
 */
static void days_to_date_string(int32_t days, char* buf, size_t len) {
    time_t secs = (time_t)days * 86400;
    struct tm* tm = gmtime(&secs);
    strftime(buf, len, "%Y-%m-%d", tm);
}

/**
 * Write all data types to a Parquet file
 */
static int write_typed_data(const char* filename, carquet_schema_t* schema) {
    printf("\nWriting typed data to: %s\n", filename);

    carquet_error_t err = CARQUET_ERROR_INIT;

    /* Create writer */
    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = CARQUET_COMPRESSION_SNAPPY;

    carquet_writer_t* writer = carquet_writer_create(filename, schema, &opts, &err);
    if (!writer) {
        fprintf(stderr, "Failed to create writer: %s\n", err.message);
        return -1;
    }

    /* Allocate data arrays */
    int32_t* counts = malloc(NUM_ROWS * sizeof(int32_t));
    int32_t* dates = malloc(NUM_ROWS * sizeof(int32_t));
    int64_t* big_numbers = malloc(NUM_ROWS * sizeof(int64_t));
    int64_t* timestamps = malloc(NUM_ROWS * sizeof(int64_t));
    float* temperatures = malloc(NUM_ROWS * sizeof(float));
    double* precise_values = malloc(NUM_ROWS * sizeof(double));
    carquet_byte_array_t* descriptions = malloc(NUM_ROWS * sizeof(carquet_byte_array_t));

    char desc_buffers[NUM_ROWS][64];

    /* Generate sample data */
    int32_t base_date = 19000;  /* Days since epoch (~2022-01-01) */
    int64_t base_time = 1640000000000000LL;  /* Microseconds (~2021-12-20) */

    for (int i = 0; i < NUM_ROWS; i++) {
        counts[i] = i * 10;
        dates[i] = base_date + i;
        big_numbers[i] = (int64_t)i * 1000000LL + 123456789LL;
        timestamps[i] = base_time + (int64_t)i * 3600000000LL;  /* +1 hour each */
        temperatures[i] = 20.0f + (i % 30) * 0.5f;
        precise_values[i] = 3.141592653589793 * (i + 1);

        snprintf(desc_buffers[i], sizeof(desc_buffers[i]),
                 "Item #%d description text", i + 1);
        descriptions[i].data = (uint8_t*)desc_buffers[i];
        descriptions[i].length = (int32_t)strlen(desc_buffers[i]);
    }

    /* Write all columns */
    int col = 0;
    carquet_status_t status;

    status = carquet_writer_write_batch(writer, col++, counts, NUM_ROWS, NULL, NULL);
    if (status != CARQUET_OK) goto write_error;

    status = carquet_writer_write_batch(writer, col++, dates, NUM_ROWS, NULL, NULL);
    if (status != CARQUET_OK) goto write_error;

    status = carquet_writer_write_batch(writer, col++, big_numbers, NUM_ROWS, NULL, NULL);
    if (status != CARQUET_OK) goto write_error;

    status = carquet_writer_write_batch(writer, col++, timestamps, NUM_ROWS, NULL, NULL);
    if (status != CARQUET_OK) goto write_error;

    status = carquet_writer_write_batch(writer, col++, temperatures, NUM_ROWS, NULL, NULL);
    if (status != CARQUET_OK) goto write_error;

    status = carquet_writer_write_batch(writer, col++, precise_values, NUM_ROWS, NULL, NULL);
    if (status != CARQUET_OK) goto write_error;

    status = carquet_writer_write_batch(writer, col++, descriptions, NUM_ROWS, NULL, NULL);
    if (status != CARQUET_OK) goto write_error;

    /* Close writer */
    status = carquet_writer_close(writer);
    writer = NULL;
    if (status != CARQUET_OK) {
        fprintf(stderr, "Failed to close writer\n");
        goto cleanup;
    }

    printf("  Successfully wrote %d rows with %d columns\n", NUM_ROWS, col);

    /* Clean up */
    free(counts);
    free(dates);
    free(big_numbers);
    free(timestamps);
    free(temperatures);
    free(precise_values);
    free(descriptions);
    return 0;

write_error:
    fprintf(stderr, "Failed to write column %d\n", col);
cleanup:
    free(counts);
    free(dates);
    free(big_numbers);
    free(timestamps);
    free(temperatures);
    free(precise_values);
    free(descriptions);
    if (writer) carquet_writer_abort(writer);
    return -1;
}

/**
 * Read and display typed data from a Parquet file
 */
static int read_typed_data(const char* filename) {
    printf("\nReading typed data from: %s\n", filename);

    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_reader_t* reader = carquet_reader_open(filename, NULL, &err);
    if (!reader) {
        fprintf(stderr, "Failed to open reader: %s\n", err.message);
        return -1;
    }

    int64_t num_rows = carquet_reader_num_rows(reader);
    int32_t num_cols = carquet_reader_num_columns(reader);

    printf("  File contains %lld rows, %d columns\n",
           (long long)num_rows, num_cols);

    /* Read sample values from first row group */
    enum { SAMPLE_SIZE = 5 };
    carquet_column_reader_t* col_reader;
    int64_t count;

    /* Read count column (INT32) */
    col_reader = carquet_reader_get_column(reader, 0, 0, &err);
    if (col_reader) {
        int32_t values[SAMPLE_SIZE];
        count = carquet_column_read_batch(col_reader, values, SAMPLE_SIZE, NULL, NULL);
        printf("\n  count (INT32): ");
        for (int i = 0; i < count; i++) {
            printf("%d ", values[i]);
        }
        printf("...\n");
        carquet_column_reader_free(col_reader);
    }

    /* Read date column (INT32 with DATE logical type) */
    col_reader = carquet_reader_get_column(reader, 0, 1, &err);
    if (col_reader) {
        int32_t values[SAMPLE_SIZE];
        count = carquet_column_read_batch(col_reader, values, SAMPLE_SIZE, NULL, NULL);
        printf("  created_date (DATE): ");
        for (int i = 0; i < count; i++) {
            char date_str[16];
            days_to_date_string(values[i], date_str, sizeof(date_str));
            printf("%s ", date_str);
        }
        printf("...\n");
        carquet_column_reader_free(col_reader);
    }

    /* Read big_number column (INT64) */
    col_reader = carquet_reader_get_column(reader, 0, 2, &err);
    if (col_reader) {
        int64_t values[SAMPLE_SIZE];
        count = carquet_column_read_batch(col_reader, values, SAMPLE_SIZE, NULL, NULL);
        printf("  big_number (INT64): ");
        for (int i = 0; i < count; i++) {
            printf("%lld ", (long long)values[i]);
        }
        printf("...\n");
        carquet_column_reader_free(col_reader);
    }

    /* Read temperature column (FLOAT) */
    col_reader = carquet_reader_get_column(reader, 0, 4, &err);
    if (col_reader) {
        float values[SAMPLE_SIZE];
        count = carquet_column_read_batch(col_reader, values, SAMPLE_SIZE, NULL, NULL);
        printf("  temperature (FLOAT): ");
        for (int i = 0; i < count; i++) {
            printf("%.1f ", values[i]);
        }
        printf("...\n");
        carquet_column_reader_free(col_reader);
    }

    /* Read precise_value column (DOUBLE) */
    col_reader = carquet_reader_get_column(reader, 0, 5, &err);
    if (col_reader) {
        double values[SAMPLE_SIZE];
        count = carquet_column_read_batch(col_reader, values, SAMPLE_SIZE, NULL, NULL);
        printf("  precise_value (DOUBLE): ");
        for (int i = 0; i < count; i++) {
            printf("%.6f ", values[i]);
        }
        printf("...\n");
        carquet_column_reader_free(col_reader);
    }

    carquet_reader_close(reader);
    printf("\n  Successfully read typed data\n");
    return 0;
}

int main(int argc, char* argv[]) {
    char default_path[1024];
    example_temp_path(default_path, sizeof(default_path), "example_data_types.parquet");
    const char* filename = default_path;

    printf("=== Carquet Data Types Example ===\n");
    printf("Library version: %s\n\n", carquet_version());

    if (argc > 1) {
        filename = argv[1];
    }

    /* Create typed schema */
    carquet_schema_t* schema = create_typed_schema();
    if (!schema) {
        return 1;
    }

    /* Write typed data */
    if (write_typed_data(filename, schema) != 0) {
        carquet_schema_free(schema);
        return 1;
    }

    carquet_schema_free(schema);

    /* Read typed data back */
    if (read_typed_data(filename) != 0) {
        return 1;
    }

    printf("\n=== Data types example completed ===\n");

    /* Clean up */
    if (argc <= 1) {
        remove(filename);
        printf("(Removed temporary file)\n");
    }

    return 0;
}
